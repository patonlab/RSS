# RSS
Radical Stability Score

This is a script to calculate Radical Stability Score for Quantum Chemistry output files. This is implemented using cclib, an hence output file supoorted by cllib can be used to calculate Radical Stability Score.

### Using the script through command lines in terminals
An Example of how to run it is shown below are:  
    ```
    python -m RSS --files *.log --output test --type 'mulliken'
    ```
### Required Modules
```
dbstep
cclib
pandas
numpy
```

## License
RSS is freely available under an [MIT](https://opensource.org/licenses/MIT) License  

## Reference
S. V., S. S.; Shree.; St. John, P. C.; Paton, R. S. A quantitative metric for organic radical stability and persistence using thermodynamic and kinetic features *Chem. Sci.*, **2021**,*12*, 13158-13166
