function [PartialPressure,F]=pCFCfromCFC(Concentration,Skip,Var,Temp,S);
    T=Temp+273.15;
    if Var==3
        a1=-82.1639;
        a2=120.152;
        a3=30.6372;
        b1=0.0293201;
        b2=-0.0351974;
        b3=0.00740056;
        F=exp(a1+a2*(100/T)+a3*log(T/100)+S.*(b1+b2*(T/100)+b3*(T/100).^2));
    elseif Var==1
        a1=-232.0411;
        a2=322.5546;
        a3=120.4956;
        a4=-1.39165;
        b1=-0.146531;
        b2=0.093621;
        b3=-0.0160693;
        F=exp(a1+a2*(100/T)+a3*log(T/100)+a4*(T/100).^2+S.*(b1+b2*(T/100)+b3*(T/100).^2));
    elseif Var==2
        a1=-220.2120;
        a2=301.8695;
        a3=114.8533;
        a4=-1.39165;
        b1=-0.147718;
        b2=0.093175;
        b3=-0.0157340;
        F=exp(a1+a2*(100/T)+a3*log(T/100)+a4*(T/100).^2+S.*(b1+b2*(T/100)+b3*(T/100).^2));
    else
        disp('Invalid input variable');
    end
    PartialPressure=Concentration./F;
end

    
    
