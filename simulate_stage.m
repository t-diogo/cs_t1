function [bit, Vresp, Vresn] = simulate_stage (nbits, nstage, np, Vp, Vn, Vref, Vcm, Cbinp, Cbinm, Cup_total, Cun_total, Gain, Voff)

    bit = zeros(np, nbits);
    Vcp = zeros(np, nbits+1);
    Vcn = zeros(np, nbits+1);

    % Calculate reference voltages -  Vref is halved each stage
    Vrp = Vcm + (Vref/2) * (1/2)^nstage;
    Vrn = Vcm - (Vref/2) * (1/2)^nstage;

    if nstage == 1
        % Bottom-Plate sampling
        Vcp(:,1) = 2*Vcm - Vp;
        Vcn(:,1) = 2*Vcm - Vn;
    else
        % Top-Plate sampling
        Vcp(:,1) = Vp;
        Vcn(:,1) = Vn;
    end
    
    for i = 1:nbits
        if nstage == 1
            bit(:,i) = (Vcn(:,i)-Vcp(:,i)) >= Voff;
            Vcp(:,i+1) = Vcp(:,i) - (Cbinp(i)/Cup_total)*Vcm + Vrn*not(bit(:,i))*(Cbinp(i)/Cup_total) + Vrp*bit(:,i)*(Cbinp(i)/Cup_total);
            Vcn(:,i+1) = Vcn(:,i) - (Cbinm(i)/Cun_total)*Vcm + Vrp*not(bit(:,i))*(Cbinm(i)/Cun_total) + Vrn*bit(:,i)*(Cbinm(i)/Cun_total);
        else
            bit(:,i) = (Vcp(:,i)-Vcn(:,i)) >= Voff;
            Vcp(:,i+1) = Vcp(:,i) - (Cbinp(i)/Cup_total)*Vcm + Vrp*not(bit(:,i))*(Cbinp(i)/Cup_total) + Vrn*bit(:,i)*(Cbinp(i)/Cup_total);
            Vcn(:,i+1) = Vcn(:,i) - (Cbinm(i)/Cun_total)*Vcm + Vrn*not(bit(:,i))*(Cbinm(i)/Cun_total) + Vrp*bit(:,i)*(Cbinm(i)/Cun_total);
        end
    end

    if nstage == 1
        Vresp = Vcm - (Vcp(:,nbits+1)-Vcn(:,nbits+1))*Gain / 2;
        Vresn = Vcm + (Vcp(:,nbits+1)-Vcn(:,nbits+1))*Gain / 2;
    else
        Vresp = Vcm + (Vcp(:,nbits+1)-Vcn(:,nbits+1))*Gain / 2;
        Vresn = Vcm - (Vcp(:,nbits+1)-Vcn(:,nbits+1))*Gain / 2;
    end
end
