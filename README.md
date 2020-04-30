# Shock_Waveform_Decomposition_Method
Code for the paper entitled "A general shock waveform and characterisation method".
This algorithm can decompose a mechanical shock measurement into several shock waveform components.
Please read the paper for more details.

## Disclaimer
This software is published for academic and non-commercial use only.

## Usage
There are two classes in this repository.

The 'acc' class contains time and data information within its properties that describes an acceleration measurement.
It is a subclass of 'timeseries' class, so all methods in 'timeseries' class can be directly applied.
It also provide some overloaded simple function, such as FFT, plot, and bandpass filter, which can be directly use to analyse mechanical shock signals.
Helps can be accessed by `help acc` in MATLAB's command window.

>A collection of methods dealing with acceleration time history signal of mechanical shocks. This class is written as a subclass of 'timeseries' class, so all 'timeseries' methods can also be used.
>
> acc Properties:  
>    Sf - Sample rate  
>    Time - Time column  
>    Data - Measured acceleration data column  
>    Length - Length of time series
>
>acc Methods:  
    acc - Constructor method to creat a acc object.  
    resample1 - Resample a time series.  
    bandpass - A bandpass filter.  
    plot - An overload plot function for acc object.  
    fft - An overload fast fourier transform function for acc object.  
    fit - Shock waveform decomposition method.  
    cwt - An overload continues wavelet transform plot.  
    dwt - An overload discrete wavelet transform plot.  
    cumtraapz - Overload numerical integration.  
    diff - Overload numerical difference.  
    extend - Extending the time series for a certain period.  

The 'acc.fit' method is for shock waveform decomposition.
For help for a specific method, please use `help acc.fit`, for example.

>SWD=fit(Acc, Name, Value) retures the object of shock waveform
 decomposition results. Optional name-value pair arguments can
 be added.
>
> Name-Value Pairs:
>    'FreSpace' - Spacing of frequenciey start points;
>    2 (default) | Positive scalar
>
>    'IniTim' - Start point of initial time.
>    0 (default) | Row vector
>
>    'TauLoc' - Start points of peak times.
>    [0.6,1,1.4] (default) | Row vector
>
>    'TwoWay' - Whether consider the situation with negative
>    peak time. '0' for 'No', and '1' for 'Yes'.
>    0 (default) | 1
>
>    'ErrTol' - The tolerance of error energy ratio.
>    0.1 (default) | Scalar between (0, 1)
>
>    'MinSW' - Minimum shock wavefrom components
>    5 (default) | Positive scalar
>
>    'XiList' - Start points of damping ratios.
>    logspace(-2,1,4) (default) | row vector
>
>    'PhiNum' - How many start points considered for phase.
>    2 (default) | Positive scalar
