#!/usr/bin/env bash
# Yield the shared GPU and memory to other people. Every INTERVAL seconds
# the guard looks for GPU work that is not ours (nvidia-smi pmon: any
# process not matching our pattern with more than SM_MIN percent of the
# streaming multiprocessors) and for low free memory. While either holds
# it pauses (SIGSTOP) our newest running job, one per cycle, down to
# MIN_RUNNING of them; once the GPU has been free of others and memory is
# back for CLEAR_N cycles it resumes (SIGCONT) the paused jobs one per
# cycle. Paused pids are kept in <log>.paused so "stop" can resume them.
#
# Usage: workflow/remote/gpu_yield.sh start <log> '<regex of our jobs>'
#        workflow/remote/gpu_yield.sh stop <log>
#        workflow/remote/gpu_yield.sh status <log>
# Environment: INTERVAL=30 SM_MIN=5 MEM_MIN=25 MEM_RESUME=35 MIN_RUNNING=1 CLEAR_N=3
# Write the pattern so that it does not match this script's own command
# line, e.g. '[c]2l_map.py' instead of 'c2l_map.py'.
set -u
cmd=${1:-}; log=${2:-}; pat=${3:-}
INTERVAL=${INTERVAL:-30}; SM_MIN=${SM_MIN:-5}; MEM_MIN=${MEM_MIN:-25}; MEM_RESUME=${MEM_RESUME:-35}
MIN_RUNNING=${MIN_RUNNING:-1}; CLEAR_N=${CLEAR_N:-3}
pidfile="$log.pid"; paused="$log.paused"

ours_pids() { pgrep -u "$USER" -f "$pat" | grep -vx "$$" | grep -vx "${BASHPID:-0}" | sort -n; }
state_of() { awk '{print $3}' "/proc/$1/stat" 2>/dev/null; }
mem_free_pct() { free | awk '/^Mem:/ {printf "%d", $7 * 100 / $2}'; }
others_on_gpu() {
  # compute processes (type C; graphics-only processes such as a browser or
  # the desktop are not GPU work) that are not ours and use more than
  # SM_MIN percent of the SMs
  local ours; ours=",$(ours_pids | paste -sd, -),"
  nvidia-smi pmon -c 1 -s u 2>/dev/null | awk -v ours="$ours" -v min="$SM_MIN" \
    '$1 !~ /^#/ && $2 ~ /^[0-9]+$/ && $3 ~ /C/ && ($4 + 0) > min && index(ours, "," $2 ",") == 0 {printf "%s(%s,sm%s) ", $NF, $2, $4}'
}

loop() {
  local clear=0 tick=0
  touch "$paused"
  while true; do
    tick=$((tick + 1))
    local now; now=$(date '+%F %T')
    local others; others=$(others_on_gpu)
    local mem; mem=$(mem_free_pct)
    local running=() stopped=()
    for p in $(ours_pids); do
      case "$(state_of "$p")" in T) stopped+=("$p");; "") ;; *) running+=("$p");; esac
    done
    if [[ -n "$others" || "$mem" -lt "$MEM_MIN" ]]; then
      clear=0
      if (( ${#running[@]} > MIN_RUNNING )); then
        local victim=${running[-1]}
        kill -STOP "$victim" && echo "$victim" >> "$paused" && \
          echo "$now PAUSE pid=$victim others=[${others:-none}] free=${mem}% running=$(( ${#running[@]} - 1 ))" >> "$log"
      elif (( tick % 10 == 0 )); then
        echo "$now BUSY others=[${others:-none}] free=${mem}% running=${#running[@]} paused=${#stopped[@]} (at MIN_RUNNING)" >> "$log"
      fi
    else
      clear=$((clear + 1))
      if (( clear >= CLEAR_N && ${#stopped[@]} > 0 && mem >= MEM_RESUME )); then
        local back=${stopped[0]}
        kill -CONT "$back" && sed -i "/^$back\$/d" "$paused" && \
          echo "$now RESUME pid=$back free=${mem}% running=$(( ${#running[@]} + 1 ))" >> "$log"
      elif (( tick % 20 == 0 )); then
        echo "$now ok free=${mem}% running=${#running[@]} paused=${#stopped[@]}" >> "$log"
      fi
    fi
    sleep "$INTERVAL"
  done
}

case "$cmd" in
  start)
    [[ -n "$log" && -n "$pat" ]] || { echo "usage: $0 start <log> '<pattern>'" >&2; exit 1; }
    if [[ -f "$pidfile" ]] && kill -0 "$(cat "$pidfile")" 2>/dev/null; then echo "already running (pid $(cat "$pidfile"))"; exit 0; fi
    setsid nohup "$0" _loop "$log" "$pat" > /dev/null 2>&1 < /dev/null &
    echo $! > "$pidfile"; echo "gpu_yield pid $! -> $log"
    ;;
  _loop) loop ;;
  stop)
    [[ -f "$pidfile" ]] && { kill "$(cat "$pidfile")" 2>/dev/null; rm -f "$pidfile"; }
    if [[ -f "$paused" ]]; then
      for p in $(cat "$paused"); do kill -CONT "$p" 2>/dev/null && echo "$(date '+%F %T') RESUME pid=$p (stop)" >> "$log"; done
      rm -f "$paused"
    fi
    echo stopped
    ;;
  status)
    echo "guard: $([[ -f "$pidfile" ]] && kill -0 "$(cat "$pidfile")" 2>/dev/null && echo running || echo not running)"
    echo "others on GPU: $(others_on_gpu)"; echo "free memory: $(mem_free_pct)%"
    [[ -f "$paused" ]] && echo "paused pids: $(paste -sd' ' "$paused")"
    tail -3 "$log" 2>/dev/null
    ;;
  *) echo "usage: $0 start|stop|status <log> ['<pattern>']" >&2; exit 1 ;;
esac
