function root = newton(f,df,p0,delta,epsilon,maxi)
%Apply Newton's method to compute the root of real-valued function f near p0
%
%Input:
%f: the object function
%df: the derivative of f
%p0: an initial guess
%delta: the tolorance error
%epsilon: the tolorance to admit a root
%maxi: the upper limit number of loops 
%
%Output:
%root: the root near p0
%

for i=1:maxi
    f0=f(p0);
    df0=df(p0);
    dp=-f0/df0;
    error=abs(dp);
    relerror=error/(abs(p0)+delta);
    if(error < delta || relerror < delta || abs(f0) < epsilon)
        break;
    end
    p0=p0+dp;
end
root=p0;
end

