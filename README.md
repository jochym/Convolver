Convolver
=========

Simple convolver tool to smooth and fit theory to experimental data. This is a simple GUI tool with quite spartan interface. The usage is also simple:

    ./convolv.py [--help|-h] [--exscale|-ex EX] [--eyscale|-ey EY] theory_data [exper_data] > output_data

You provide the theoretical data to convolve and experimental data to compare. The code will output your convolved data to standard output. You can use Gauss or Lorentzian shape of the convolution kernel.  

