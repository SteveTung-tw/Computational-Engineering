function out = SD_fft(fkk_original,pow_NN) %fkk_original is f(t) and the number of sampling point is 2 to the pow_NN th power
NN = 2^(pow_NN);
f0F0 = zeros(NN,pow_NN);
xx = linspace(-1*pi,pi,NN+1);


%rearrange index of fk
fks_ind = [0:NN-1] + NN;
fks_ind = dec2bin(fks_ind);
fks_ind=str2num(fliplr(num2str(fks_ind)));
fks_ind = fks_ind - 1;
fks_ind = num2str(fks_ind);
fks_ind = bin2dec(fks_ind);
fks_ind = fks_ind ./2;

fks = [];
for ii = 1:NN
    fks(ii) = fkk_original(fks_ind(ii)+1);
end


wjn = zeros(NN,pow_NN);
for ii = 1:pow_NN 
    for jj = 1:(NN/(2^ii))
        for kk = 1:(2^ii)
            wjn((jj-1)*(2^ii)+kk,ii) = exp(-1*i*2*pi*(kk-1)/(2^ii));
        end
    end
end

%%

ee_new_all = [];
oo_new_all = [];
eeoo_old = fks';

for ii = 1:pow_NN
    
    ee_new = eeoo_old;
    oo_new = eeoo_old;
    
    for jj = 1:NN/(2^(ii-1))
        
        if mod(jj,2) == 1
            copied_ee = eeoo_old((2^(ii-1))*(jj-1)+1:(2^(ii-1))*jj);
        end
        
        if mod(jj,2) == 0
            ee_new((2^(ii-1))*(jj-1)+1:(2^(ii-1))*jj) = copied_ee;
        end
        
    end
    
    for jj = NN/(2^(ii-1)):-1:1
        
        if mod(jj,2) == 0
            copied_oo = eeoo_old((2^(ii-1))*(jj-1)+1:(2^(ii-1))*jj);
        end
        
        if mod(jj,2) == 1
            oo_new((2^(ii-1))*(jj-1)+1:(2^(ii-1))*jj) = copied_oo;
        end
        
    end
    ee_new_all = [ee_new_all,ee_new];
    oo_new_all = [oo_new_all,oo_new];
    eeoo_new = (ee_new .* 1) + (oo_new .* wjn(:,ii));
    f0F0(:,ii) = eeoo_new;
    eeoo_old = eeoo_new;
end


F0 = f0F0(:,end);
F0 = F0';
out = F0;