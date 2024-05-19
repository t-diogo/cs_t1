function [dout] = simulate_pipeline (num_stages, nbits, Vref, Vcm, deltaVin, np, Cbinp, Cbinm, Cup_total, Cun_total, Gain, Voff)

    dout_stages = cell(1,num_stages);

    % Input signal conversion - differential to single ended
    Vp = Vcm + deltaVin/2;
    Vn = Vcm - deltaVin/2;

    % Pipeline simulation
    for i = 1:num_stages
        [dout_stages{i}, Vp, Vn] = simulate_stage(nbits(i), i, np, Vp, Vn, Vref, Vcm, Cbinp{i}, Cbinm{i}, Cup_total{i}, Cun_total{i}, Gain{i}, Voff{i});
    end

    dout1 = binaryVectorToDecimal(dout_stages{1});
    dout2 = binaryVectorToDecimal(dout_stages{2});
    dout3 = binaryVectorToDecimal(dout_stages{3});
    
    % digital correction
    dout = dout1*(2^8) + (dout2 - 2^(nbits(2)-2))*(2^5)  + (dout3 - 2^(nbits(3)-2));

end