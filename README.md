<H1>Matlab codes for potential theoretic design of formulas for approximation of functions</H1>

These are Matlab R2015a (8.5.0.197613) codes for potential theoretic design of formulas for approximation of functions in weighted Hardy spaces. These codes were used to yield the results of the numerical experiments in the article

"Potential theoretic approach to design of highly accurate formulas for function approximation in weighted Hardy spaces" (in preparation),

which was written by my collaborators and me. See this article for the technical details of the method.

The roles of the source files and their relationships are explained below.

[Main program for the design of the proposed formulas]

<ul>
  <li> MAIN_opt_sample_wrt_N.m &middot;&middot;&middot; A program to generate the sampling points of the proposed formulas. The subroutines used in this program are listed below.
  <ul>
    <li> SUB_alpha.m &middot;&middot;&middot; A subroutine to obtain the approximate value of &alpha;<sup>*</sup><sub>N</sub> with the Newton method. </li>
    <li> SUB_gen_opt_sample.m &middot;&middot;&middot; A subroutine to generate the sampling points with the Fourier transform. The following auxiliary programs are used. The first two are used to define some mathematical functions and the last one executes fractional FFT.</li>
    <ul>
      <li> SC.m </li>
      <li> FLT.m </li>
      <li> FFFT.m </li>
    </ul>
  </ul>
  </li>

</ul>

[Main program for the experiments of approximating functions by the proposed formulas]

