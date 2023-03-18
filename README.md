# model_slip_rate_GPS

MATLAB script to model the Slip Rate and Locking Depth of a fault from GPS data

Matrix Inversion uses the form A*x = b

The matrix A contains the equation (1/pi)*atan(x/D) as first column and 1 as second column
Note the matrix x contains the slip rate value and constant which we need to find out for different locking Depth values.
Matrix b contains the velocity values.
