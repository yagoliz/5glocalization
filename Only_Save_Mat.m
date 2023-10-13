clear all
close all

pos_UE_vector = [1.80,6.07;1.00,6.07;1.80,9.14;3.28,2.97;0.58,2.03;3.28,7.68];

pos_gNB_1 = NaN(2,6);
pos_gNB_2 = NaN(2,6);
pos_gNB_3 = NaN(2,6);

ToA_All = NaN(4,5,100);

compensation_1 = 0;
compensation_2 = 0;
compensation_3 = 0;

for reference_pos =  5%:6

    % gNBs position
    XR1=3.85;
    YR1=12.83;
    XR2=0;
    YR2=12.83;
    XR3=0;
    YR3=0;
    XR4=3.85;
    YR4=0;

    bw=80;
    id_meas=1;
    disp(strcat('experiments/exp',int2str(reference_pos-1),'_',int2str(bw),'mhz_',int2str(id_meas),'.txt'))

    fid = fopen(strcat('experiments/exp',int2str(reference_pos-1),'_',int2str(bw),'mhz_',int2str(id_meas),'.txt'))


    tline = fgetl(fid);
    read_measure = 1;
    toa = [];
    which_gNB = [];
    while ischar(tline)
        toa_index = strfind(tline,'DL PRS ToA ==>');
        end_toa_index =strfind(tline,', peak');
        toa(read_measure)=str2num(tline((toa_index+15):(end_toa_index-16)));

        gNB_index = strfind(tline,'gNB');
        which_gNB(read_measure)=str2num(tline(gNB_index+4));
        
        sfn_index = strfind(tline,'sfn');
        V=[3:10];
        end_sfn=strfind(tline(sfn_index+V),']');
        which_sfn(read_measure)=str2num(tline((sfn_index+3):(sfn_index+V(end_sfn)-1)));
        tline = fgetl(fid);
        tline = fgetl(fid);
        
        read_measure=read_measure+1;
    end
    
    % matrix ToA
    tot_measure=floor(length(toa)/4)*4;
    ToA_Matrix=(reshape(toa(1:tot_measure),4,[]));
    gNB_Matrix=reshape(which_gNB(1:tot_measure),4,[]);
    
  %  ToA_Comp = [];
    
    pos_gNB_1(:,reference_pos) = [XR1,YR1];
    pos_gNB_2(:,reference_pos) = [XR2,YR2];
    pos_gNB_3(:,reference_pos) = [XR3,YR3];
    pos_gNB_4(:,reference_pos) = [XR4,YR4];
    
%     pos_UE_vector
    
    %ToA_All(:,reference_pos,1:size(ToA_Matrix,2)) = ToA_Matrix;
    
    % SNR, RSRP, peak
    for gg=1:4
        clearvars TT
        TT(:)=ToA_Matrix(gg,:);
        hold on;
        median(TT)  
        histogram(TT);
    end
end

sqrt(sum((pos_gNB_1-pos_UE_vector').^2))

