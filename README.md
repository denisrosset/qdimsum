# QDimSum

Symmetric SDP relaxations for qudits systems

The code requires YALMIP and a semidefinite solver of your choice (for example Mosek or SeDuMi).

See [RAC22.m](https://github.com/denisrosset/qdimsum/blob/master/RAC22.m) and [TestRAC22.m](https://github.com/denisrosset/qdimsum/blob/master/TestRAC22.m) for a simple example.


## Limitations

- Code should work on complex moment matrices, but this has not been tested. (Please contact us with a use case!).
- Only decomposition of representations when all irreps are real (complex and quaternionic irreps work in progress).
