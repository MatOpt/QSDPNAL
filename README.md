## QSDPNAL version 0.1 -- a MATLAB software for convex quadratic semidefinite programming

### [Xudong Li](http://lixudong.info), [Defeng Sun](https://www.polyu.edu.hk/ama/profile/dfsun/), [Kim-Chuan Toh](https://blog.nus.edu.sg/mattohkc/)

This software is designed to solve convex quadratic SDP of the form:
    
    min 1/2<X,QX> + <C,X> subject to A(X)=b, B(X)>=d, L <= X <= U, X is PSD

### Citation

- **Xudong Li, Defeng Sun, and Kim-Chuan Toh**, QSDPNAL: A two-phase augmented Lagrangian method for convex quadratic semidefinite programming, Mathematical Programming Computation, 10 (2018), pp. 703--743. 

-------------

- **Important note**:  
  - The software is still under development. Thus it will invariably be buggy. We would appreciate your feedback and bugs’ report.
  - This is a research software. It is not intended nor designed to be a general purpose software at the moment.
  
-----

### Copyright
This version of QSDPNAL is distributed under the **3-Clause BSD license**. For commercial applications that may be incompatible with this license, please contact the authors to discuss alternatives.

-----------

### Installation 

Welcome to QSDPNAL! Users can simply follow the steps below to install QSDPNAL within Matlab:

  - Firstly, clone the package via
    ```github
    git clone https://github.com/MatOpt/QSDPNAL.git
    ```
  - Run Matlab in the directory QSDPNAL 
  - In Matlab command window, type:  
    ```matlab
    >> startup
    >> qsdpnaldemo
    ```
  - <font color=blue>By now, QSDPNAL is ready for you to use</font>.
  - If you have issues with MEX files on macOS (“ “*.mexmaci64” cannot be opened because the developer cannot be verified. macOS cannot verify that this app is free from malware” or “Code signature not valid for use in process using Library Validation: library load disallowed by system policy”), open a Terminal, cd to the QSDPNAL directory and type:
  
    `find . -name "*.mexmaci64" -exec xattr -d com.apple.quarantine {} \;`

  - The **User's guide** is included in the [Appendix of our paper](https://link.springer.com/article/10.1007/s12532-018-0137-6).

---------------
  - Codes in the folder QSDPNAL_main\NAL are for Algorithm QSDPNAL (QSDPNAL-Phase II with QSDPNAL-Phase I to generate initial point)

  - Codes in the folder QSDPNAL_main\SCB are for Algorithm QSDPNAL-Phase I (SCB_ADMM: an SCB based inexact semi-proximal ADMM)

  - Codes in the folder QSDPNAL_main\submain\SNCG are used for solving subproblems (i.e., Algorithms SSNCG and ABCD)
