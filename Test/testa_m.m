h=0.1;
m=1;
R1=1.25;
R2=1.25;
R3=1.25;
R4=1.25;
R5=1.25;
l=5*pi/3;
k=1;
mode=1;
d1=sqrt(3)/2;
Out=[];
format long;
for hs=[h,h/4,h/4.^2,h/4.^3]
    a_m1=asymp_a_m(m,k,l,R1,hs,d1,mode);
    a_m2=asymp_a_m(m,k,l,R2,hs,d1,mode);
    a_m3=asymp_a_m(m,k,l,R3,hs,d1,mode);
    a_m4=asymp_a_m(m,k,l,R2,hs,d1,mode);
    a_m5=asymp_a_m(m,k,l,R3,hs,d1,mode);
    a_m6=TE_RHS5(m,0,l,k,R1,R2,R3,R4,R5,hs,hs,hs,hs,hs,mode,d1);
    %if(mode==1)
    %    a_m2=a_m2/1i;
    %end
    delta=[a_m1,a_m2,a_m3,a_m4,a_m5]-a_m6.';
    Out=[Out;[a_m1,a_m2,a_m3,a_m4,a_m5],a_m6.',delta];
end
display(Out);


Out2=[];
format long;
for hs=[h,h/4,h/4.^2,h/4.^3]
    a_m3=TE_RHS(m,5,1,k,R1,hs,mode,d1);
    delta=a_m3(2:11);
    Out2=[Out2;((delta(1:end)).'*delta(1:end)).^(1/2)];
end
display(Out2);