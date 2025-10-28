% MOD_LA_SML

% % Convert from to mol/(km^2 hr) to molecules/(cm^3 s) using the formula:
% (6.022e23 molecules/mol) * (1/BLH(m)*1000m/1km) * (1km^3/1e15cm^3) * (1hr/3600s)
conv_fac = 1.67278e8 ./ BLH;

i=i+1;
Rnames{i} = 'Emis = APINENE';
k(:,i) = EAPI .* conv_fac;  % Read from Met in molec/cm^3/s
fAPINE(i) = 1;  % Ensure 100% of emitted mass goes to α-pinene

i=i+1;
Rnames{i} = 'Emis = BPINENE';
k(:,i) = EBPI .* conv_fac;
fBPINE(i) = 1;

i=i+1;
Rnames{i} = 'Emis = LIMONENE';
k(:,i) = ELIM .* conv_fac;
fDLIMO(i) = 1;

i=i+1;
Rnames{i} = 'Emis = ISOPRENE';
k(:,i) = EISP .* conv_fac;
fISOP(i) = 1;
