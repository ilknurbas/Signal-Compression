% İlknur Baş
clc; close all; clear; 
%% (1)
% name of the compressed image is 'BITSTREAM_FILE.bin'

A = imread('Image_4.png'); % 1,2,3 ,4
A = double(A);

figure(1)
histogram(A(:));
xlabel('Gray values');
ylabel('Frequencies');
title('The plot of the histogram for the image 4');

figure(2)
histogram(A(:), 'Normalization', 'probability');
xlabel('Gray values');
ylabel('Frequency Probabilities');
title('The plot of the histogram for the image 4');

%% (2) 
B = A(51:306,51:306);
size(B); % 256   256
imwrite(uint8(B),'MyImage_4.png');
sd = dir('MyImage_4.png');
filesize = sd.bytes; % 39.992;
size(B(:)); % 65.536           1

%%
figure(3)
histogram(B(:),'Normalization', 'count');
xlabel('Gray values');
ylabel('Frequencies');
title('The plot of the histogram for the image B');

figure(4)
histogram(B(:),'Normalization', 'probability');
xlabel('Gray values');
ylabel('Frequency Probabilities');
title('The plot of the histogram for the image B');
%% 
hist_prob_b = histogram(B(:), 'Normalization', 'probability');
probabilities_found = hist_prob_b.Values(hist_prob_b.Values ~= 0);

entropy = -sum(probabilities_found.*log2(probabilities_found));
disp(['Empirical entropy for the image B: ', num2str(entropy)]);

%% (3)
E = zeros(size(B));
[row,column] = size(B);
for r = 1:row
    for c = 1:column
        if r==1 && c==1
            y_hat = 128;
        else
            n = nan;
            w = nan;
            nw = nan;

            if r>1
                n = B(r-1,c);
            end 

            if c>1
                w = B(r,c-1);
            end 

            if r>1 && c>1
                nw = B(r-1,c-1);
            end 
            
            inside = [n w n+w-nw];
            inside = inside(~isnan(inside));
            y_hat = round(median(inside));
            
        end

        E(r,c) = B(r,c)-y_hat;
    end
end

figure(5);
imagesc(E); % E(i, j) are all integers belonging to the set {−255, . . . , 255}.

% shifting (not much difference)
E_shift = (E + 255) / 510; % Scale elements to the range of 0 to 1
figure(6);
imagesc(E_shift);  % Use imagesc for visualizing a matrix with color mapping

% concatenate columnwise
e = E(:)';
size(e); % 1       65536

%% (6)
L = 50:50:2000; % lenght 40
CL = zeros(1, length(L));
for i=1:length(L)
    CL(i) = Segment_L(e,L(i));
end

% Find the block size that produces the smallest codelength and 
% use it to finally encode the vector evec.
CL;
[min_codelenght, ind] = min(CL);
min_codelenght;
ind;
min_block_size = L(ind);
disp(['The block size that produces the smallest codelength:', num2str(min_block_size)]);

% Compute bpp CL using CLi by dividing each value in CLi with 
% the number of pixels existing in the matrix B. Hence, bpp CL contains 
% the number of bits per pixel needed to compress the matrix B.
bpp_CL = CL/length(B(:));
bpp_CL;

figure(7);
plot(L,bpp_CL);
xlabel('L values - block size')
ylabel('bits per pixel')
title('the plot of bpp CLi versus L')

% (4)
% encoding
disp('Encoding starts...');
bitstream = GR_encode(e,min_block_size);
fid = fopen('BITSTREAM_FILE.bin', 'w');
fwrite(fid, bitstream, 'ubit1');
fclose(fid);
disp("Encoding is done.");

% (5)
% decoding
disp("Decoding starts...");
fid = fopen('BITSTREAM_FILE.bin', 'r');
bitstream = fread(fid, inf, 'ubit1');
fclose(fid);
decoded_e = decode_block(bitstream, min_block_size);
disp("Decoding is done.");

not_same_indices = find(e ~= decoded_e);
num_not_same = length(not_same_indices);
size(decoded_e);
disp(['Number of indices that are not the same (expected 0): ' , num2str(num_not_same)]);
disp(['Error in reconstruction (expected 0): ' num2str(sum((e-decoded_e)))]);

decoded_B = reshape(decoded_e, 256, 256);
figure(8)
imagesc(decoded_B)
title('decoded B matrix')

% Functions
function [decoded_e] = decode_block(bitstream,min_block_size)
    
    array_length = 65536;
    decoded_e = zeros(1, array_length);
    bitstram_start_counter = 1; 
    end_block = bitstram_start_counter;
    min_block_size;
    for i = 1:min_block_size:array_length
         
        if i+min_block_size > array_length
            ends = array_length;
        else
            ends = i + min_block_size - 1;
        end
    
        [block_e, end_block] = decode_single(bitstream, end_block, ends-i+1);
         
        decoded_e(i:ends) = block_e;
    end

end

function [block_e, end_block] = decode_single(bitstream, counter, block_lenght)

    % first 4 has the p value
    starts = counter;
    ends = starts+4-1;
    (bitstream(starts:ends));
    p = bi2de((bitstream(starts:ends).'),'left-msb'); % ,'left-msb'
    starts = ends+1;
    
    % signs
    signs = zeros(block_lenght, 1);
    quotients = zeros(block_lenght, 1);
    remainders = zeros(block_lenght, 1);

    if p == 0 
        for i = 1:block_lenght

            % quotient
            % find 1
            end_1 = starts;
            
            while (bitstream(end_1) == 0)
                bitstream(end_1);
                end_1;
                end_1 = end_1 + 1;
            end
            end_1;
            unary_bin = bitstream(starts:end_1);
            quotient = length(unary_bin) - 1;
    
            starts = end_1 + 1;
    
            % sign
            % sign = bitstream(starts:starts);
            sign = bitstream(starts);
            starts = starts+1;
    
            % save
            signs(i) = sign;
            quotients(i) = quotient;
            remainders(i) = 0;
    
        end
    
    else  
        for i = 1:block_lenght
            
            % remainder bit --p
            ends = starts + p-1;
            remainder_bits = bitstream(starts:ends);
            remainder = bi2de(remainder_bits.','left-msb');
       
            starts = ends +1;
    
            % quotient
            % find 1
            end_1 = starts;
            
            while (bitstream(end_1) == 0)
                bitstream(end_1);
                end_1;
                end_1 = end_1 + 1;
            end
            end_1;
            unary_bin = bitstream(starts:end_1);
            quotient = length(unary_bin) - 1;
    
            starts = end_1 + 1;
    
            % sign
            % sign = bitstream(starts:starts);
            sign = bitstream(starts);
            starts = starts+1;
    
            % save
            signs(i) = sign;
            quotients(i) = quotient;
            remainders(i) = remainder;
    
        end

    end
    
    end_block = starts;
    % 0 to the bitstream, signalling that x is positive
    signs(signs == 1) = -1;
    signs(signs == 0) = 1;

    block_e = signs.*(quotients*(2^p)+remainders);
 
    
end

% this function encodes "values"
function [bitstream] = GR_encode(values,min_block_size)
    array_length = length(values); % 65536
    % find p values according to the min_block_size
    [p_star, total_cl] = find_optimal_p(values,min_block_size);
    p_star; % lenght 32 and 2
    total_cl;
    % bitstream = []; % initialize
    bitstream = zeros(total_cl,1);
    total_cl; % 329129

    counter = 1;
    start = 1;
    for i = 1:min_block_size:array_length
        if i+min_block_size > array_length
            block_values = values(i:end);
        else
            block_values = values(i:i+min_block_size-1);
        end
    
        block_compressed = compress_block(block_values, p_star(counter));
        counter = counter + 1;
    
        % append to bitstream
        ends = start+length(block_compressed)-1;
        bitstream(start:ends) = block_compressed;
        start = ends+1;
    
    end

end


function [block_compressed_bitstream] = compress_block(blc_values,pstar_single)
    % note : codelenght = sum(p + floor(values/2^p + 1) + 1) + 4;
    codelength = GR_estimation(blc_values,pstar_single);
    array_length = length(blc_values);
    block_compressed_bitstream = zeros(1, codelength);
     
    if pstar_single < 15
        get_bits = bitget(pstar_single, 4:-1:1) ;
    else
        disp('*****pstar_single > 15')
        %get_bits = bitget(pstar_single, 4:-1:1) ;
    end
    
    %if pstar_single == 0 
    %    disp('*****p == 0 ')
    %end

    starts = 1;
    ends =  starts+4-1;
    % You should add in the end a 4 to CL(X,p) 
    % to account for the need to transmit to the decoder the value of p that it needs to use.
    block_compressed_bitstream(starts:ends) = get_bits;
    block_compressed_bitstream(starts:ends) ;

    starts = ends+1;
    
    % Finally since the number x is non-zero, 
    % we transmit its sign bit (in this case we append a 0 to the 
    % bitstream, signalling that x is positive).
    signs = zeros(size(blc_values));
    signs(blc_values < 0) = 1;
    
    if pstar_single == 0  % gave error for images 2-3
        quotient = (abs(blc_values))/(2^pstar_single); % no double
        
        for i = 1:array_length
            single_quotient = quotient(i);
        
            %quotient_unary = [repmat('0', 1, single_quotient), '1'];
            quotient_unary = [zeros(1, single_quotient) 1];
    
            sign_bits = signs(i);
            
            % starts sign(1) pstar_single  single_quotient unary(1)
            ends = starts+1+single_quotient+1-1;
            starts; % 5
            ends; % 20
            (quotient_unary); % 12
            sign_bits; % 1

            combined_bits = [quotient_unary sign_bits];
            %a = num2str(combined_bits);

            quotient_unary;
            block_compressed_bitstream(starts:ends) = combined_bits;    
            block_compressed_bitstream(starts:ends);
    
             
            starts = ends+1;
        end
    else 
        remainder = mod(abs(blc_values), 2^pstar_single);
        remainder_bits = de2bi(remainder, pstar_single,'left-msb');
        quotient = (abs(blc_values)-remainder)/(2^pstar_single); % no double
        
        for i = 1:array_length
            single_quotient = quotient(i);
        
            %quotient_unary = [repmat('0', 1, single_quotient), '1'];
            quotient_unary = [zeros(1, single_quotient) 1];
    
            single_remainder_bits = remainder_bits(i, :);
            sign_bits = signs(i);
            
            % starts sign(1) pstar_single  single_quotient unary(1)
            ends = starts+1+pstar_single+single_quotient+1-1;
            starts; % 5
            ends; % 20
            (quotient_unary); % 12
            sign_bits; % 1
            single_remainder_bits; % 2
            length([single_remainder_bits quotient_unary sign_bits]); % 15
            combined_bits = [single_remainder_bits quotient_unary sign_bits];
            %a = num2str(combined_bits);
            single_remainder_bits;
            quotient_unary;
            block_compressed_bitstream(starts:ends) = combined_bits;    
            block_compressed_bitstream(starts:ends);
    
             
            starts = ends+1;
        end

    end


end


% L denotes a single block_size
function [codelenght_e_vec] = Segment_L(e_vec,L)
    [~, total_cl] = find_optimal_p(e_vec,L);
    codelenght_e_vec = total_cl;
    
end

% p_star returns best p value for each block
% total_cl returns total codelenght for input "values" for each block having 
% specified "block_size"
% in our case "values" variable is e 
function [p_star, total_cl] = find_optimal_p(values,block_size)

    p_star = zeros(1,ceil(length(values)/block_size));
    no_block = 1;
    total_cl = 0;
    for i=1:block_size:length(values)
        % find the values in each block
        % if its not divisible by block_size
        if i+block_size >length(values)
            values_block = values(i:end);
        else
            %display(i)
            %display(i+block_size-1)
            size(values);
            values_block = values(i:i+block_size-1);
        end
        % store codelenght values for each p for a single block
        cl = zeros(1,9);
        for p_values = 0:8
            codelenght = GR_estimation(values_block,p_values);
            cl(p_values+1) = codelenght;
        end

        [m, index] = min(cl);
        p_star(no_block) = index-1;
        no_block = no_block + 1;
        total_cl = total_cl + m;

    end

end

function [codelenght] = GR_estimation(values,p)
    values = abs(values); % abs is needed for sum
    delta = 1; % always
    codelenght = sum(p + floor(values/2^p + 1) + delta) + 4;
end





