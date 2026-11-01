# pnatlas workflow. Rules are grouped by stage under workflow/rules/;
# docs/workflow.md will describe the order, inputs and outputs of each step.
configfile: "config/config.yaml"
configfile: "config/paths.yaml"


rule all:
    default_target: True
    input: []
