#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool

baseCommand: [./base_module.py]

inputs:

  command:
    type: File
    default: 
      class: File
      path: ./base_module.py

  first_number:
    type: int
    inputBinding:
      position: 1
      prefix: --firstnum
    label: "first-number"
    doc: "first number to add"

  second_number:
    type: int
    inputBinding:
      position: 2
      prefix: --secondnum
    label: "second-number"
    doc: "second number to add"

stdout: output.txt
# stdout: $(inputs.first_number.nameroot)_plus_$(inputs.second_number.nameroot).txt

outputs:

  output_file:
    type: File
    outputBinding:
      glob: output.txt # $(inputs.first_number.nameroot)_plus_$(inputs.second_number.nameroot).txt
    label: ""
    doc: "File containing output of the base program"
    

