Background
=============
JAMA is a basic linear algebra package for Java. It provides user-level classes for constructing and manipulating real, dense matrices. It is meant to provide sufficient functionality for routine problems, packaged in a way that is natural and understandable to non-experts. It is intended to serve as the standard matrix class for Java, and will be proposed as such to the [Java Grande Forum][1] and then to [Sun][2]. A straightforward public-domain reference implementation has been developed by the [MathWorks][3] and [NIST][4] as a strawman for such a class. We are releasing this version in order to obtain public comment. There is no guarantee that future versions of JAMA will be compatible with this one.

A sibling matrix package, [Jampack][5], has also been developed at NIST and the University of Maryland. The two packages arose from the need to evaluate alternate designs for the implementation of matrices in Java. JAMA is based on a single matrix class within a strictly object-oriented framework. Jampack uses a more open approach that lends itself to extension by the user. As it turns out, for the casual user the packages differ principally in the syntax of the matrix operations. We hope you will take the time to look at Jampack along with JAMA. There is much to be learned from both packages.

Capabilities
=============
JAMA is comprised of six Java classes: Matrix, CholeskyDecomposition, LUDecomposition, QRDecomposition, SingularValueDecomposition and EigenvalueDecomposition.
The Matrix class provides the fundamental operations of numerical linear algebra. Various constructors create Matrices from two dimensional arrays of double precision floating point numbers. Various gets and sets provide access to submatrices and matrix elements. The basic arithmetic operations include matrix addition and multiplication, matrix norms and selected element-by-element array operations. A convenient matrix print method is also included.

There are methods to apply elementary functions (e.g. Math::sin, Math::exp, Math::log) on a single matrix elementwise, as well as to apply binary function on two matrices elementwise (both called transform). They can easily be used with lambda expressions in Java 8.

Five fundamental matrix decompositions, which consist of pairs or triples of matrices, permutation vectors, and the like, produce results in five decomposition classes. These decompositions are accessed by the Matrix class to compute solutions of simultaneous linear equations, determinants, inverses and other matrix functions. The five decompositions are:
* Cholesky Decomposition of symmetric, positive definite matrices
* LU Decomposition (Gaussian elimination) of rectangular matrices
* QR Decomposition of rectangular matrices
* Eigenvalue Decomposition of both symmetric and nonsymmetric square matrices
* Singular Value Decomposition of rectangular matrices

The current JAMA deals only with real matrices. We expect that future versions will also address complex matrices. This has been deferred since crucial design decisions cannot be made until certain issues regarding the implementation of complex in the Java language are resolved.
The design of JAMA represents a compromise between the need for pure and elegant object-oriented design and the need to enable high performance implementations.

Reference Implementation
=============
The implementation of JAMA downloadable from this site is meant to be a reference implementation only. As such, it is pedagogical in nature. The algorithms employed are similar to those of the classic Wilkinson and Reinsch Handbook, i.e. the same algorithms used in [EISPACK][6], [LINPACK][7] and [MATLAB][3]. Matrices are stored internally as native Java arrays (i.e., double[][]). The coding style is straightforward and readable. While the reference implementation itself should provide reasonable execution speed for small to moderate size applications, we fully expect software vendors and Java VMs to provide versions which are optimized for particular environments.

Not Covered
=============
JAMA is by no means a complete linear algebra environment. For example, there are no provisions for matrices with particular structure (e.g., banded, sparse) or for more specialized decompositions (e.g. Shur, generalized eigenvalue). Complex matrices are not included. It is not our intention to ignore these important problems. We expect that some of these (e.g. complex) will be addressed in future versions. It is our intent that the design of JAMA not preclude extension to some of these additional areas.
Finally, JAMA is not a general-purpose array class. Instead, it focuses on the principle mathematical functionality required to do numerical linear algebra.


Authors
=============
JAMA's initial design, as well as this reference implementation, was developed by
* Joe Hicklin
* [Cleve Moler](https://www.mathworks.com/company/aboutus/founders/clevemoler.html)
* Peter Webb 

...from [The MathWorks][3]

* [Ronald F. Boisvert](https://www.nist.gov/people/ronald-f-boisvert)
* [Bruce Miller](https://www.nist.gov/people/bruce-r-miller)
* [Roldan Pozo](http://math.nist.gov/~RPozo/)
* [Karin Remington](http://math.nist.gov/~KRemington/)

...from [NIST][4]

Copyright Notice
=============
This software is a cooperative product of The MathWorks and the National Institute of Standards and Technology (NIST) which has been released to the public domain. Neither The MathWorks nor NIST assumes any responsibility whatsoever for its use by other parties, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. 


Related Links &amp; Libraries
=============
* [Java for Computational Science and Engineering](http://www.npac.syr.edu/projects/javaforcse/)
* [Java Numerics Working Group](http://math.nist.gov/javanumerics/)
* [NIST Mathematical and Computational Sciences Division](http://math.nist.gov/mcsd/)

As Jama is in the public domain, other developers are free to adopt and adapt this code to other styles of programming or to extend or modernize the API.
You might find one of these libraries to be more suitable to your purposes.
Make note, however, that NIST makes no endorsement of these projects.  We are currently aware of the following ports of Jama:

* [de.mukis.jama](https://muuki88.github.com/jama-osgi/) is hosted at Apache's Maven, with unit tests.

[1]: http://www.npac.syr.edu/javagrande/
[2]: http://java.sun.com
[3]: http://www.mathworks.com
[4]: http://www.nist.gov
[5]: ftp://math.nist.gov/pub/Jampack/Jampack/AboutJampack.html
[6]: http://www.netlib.org/eispack/
[7]: http://www.netlib.org/linpack/
