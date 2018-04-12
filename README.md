Biquad Filter Class
===================

This is a simple biquad filter class that enables live digital filtering on real-time devices and microcontrollers or signal processing on all other computer devices.

Also check out [tomlankhorst/control](https://github.com/tomlankhorst/control)!

The file main.cpp contains an application example.

![Step response of two digital Butterworth low-pass filters](https://tomlankhorst.nl/wp-content/uploads/2016/09/stepresponse.png)

Please refer to [my blog post](https://tomlankhorst.nl/filter-controller-cpp-implementation/) for more details on the subject and application.

Generating C++ code from a MATLAB transfer-function
--------

The following MATLAB function converts a SOS matrix to C++ code:
```MATLAB
function [ ] = tf2cppbq( sos )
%TF2CPPBQ( sos ) Transfer-function to C++ code that initializes BiQuads and BiQuad chain
% Input: matrix of second-order-sections (use tf2sos(H) for example).

fprintf('\n');

i = 0;
for s = sos.'
    i = i + 1;
    fprintf('BiQuad bq%d( %.5e, %.5e, %.5e, %.5e, %.5e );\n', i, s(1), s(2), s(3), s(5), s(6));
end

fprintf('\nBiQuadChain bqc;\n');
fprintf('bqc');
i = 0;
for s = sos.'
    i = i + 1;
    fprintf('.add( &bq%d )', i);
end

fprintf(';\n');

end
```

Use it as follows:

```MATLAB
[b,a] = butter(4,0.2,'low');
sos   = tf2sos(b,a);
tf2cppbq( sos );
```

Output:
```C++
BiQuad bq1( 4.82434e-03, 9.64869e-03, 4.82434e-03, -1.04860e+00, 2.96140e-01 );
BiQuad bq2( 1.00000e+00, 2.00000e+00, 1.00000e+00, -1.32091e+00, 6.32739e-01 );

BiQuadChain bqc;
bqc.add( &bq1 ).add( &bq2 );
```


License
-------
MIT
