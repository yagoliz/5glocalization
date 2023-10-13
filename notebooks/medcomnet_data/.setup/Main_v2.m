clear all
close all

compensation_1 = [0];
compensation_2 = [0];
compensation_3 = [0];

for iii = 1:1
    
    pos_UE_vector = [3.3,1.6 ; 7.03,2.83 ; 7,3 ; 5.6,3 ; 4.2,3 ; ...
        3,3 ; 1.7,3 ; 1.7,1.9 ; 3.4,1.9 ; 4.5,1.9; ...
        4.4,3.55 ; 3.4,3.55 ; 2.2,3.5 ; 3,3 ; 5.6,3 ; ...
        7,3 ; 8.8,3.3];
    
    sliding_windows = 10;
        
    pos_err = [];
    
    doa12_error = [];
    doa13_error = [];
    doa23_error = [];
    tdoa12_error = [];
    tdoa13_error = [];
    tdoa23_error = [];
    toa1_error = [];
    toa2_error = [];
    toa3_error = [];
    
    for reference_pos =  1:17
        
        % gNBs position
        if reference_pos > 2
            XR1=11.52;
            YR1=3;
            XR2=0;
            YR2=5.1;
            XR3=0;
            YR3=0;
        else
            XR1=11.35;
            YR1=0;
            XR2=1.3;
            YR2=4.8;
            XR3=0;
            YR3=0;
        end
        if reference_pos==16
            XR1=15.6;
            YR1=3;
        end
        
        % Open file
        if reference_pos < 14
            fid = fopen(strcat('exp',int2str(reference_pos),'/prs_rx_',int2str(reference_pos)));
        else
            fid = fopen(strcat('exp',int2str(reference_pos),'/rep4/prs_rx_',int2str(reference_pos)));
        end
        
        pos_UE=pos_UE_vector(reference_pos,:);
        
        tline = fgetl(fid);
        read_measure = 1;
        toa = [];
        which_gNB = [];
        while ischar(tline)
            toa_index = strfind(tline,'dl_toa: ');
            end_toa_index =strfind(tline,', peak');
            toa(read_measure)=str2num(tline((toa_index+8):(end_toa_index-1)));
            which_gNB(read_measure)=str2num(tline(9));
            tline = fgetl(fid);
            read_measure=read_measure+1;
        end
        
        % matrix ToA
        tot_mesure=floor(length(toa)/3)*3;
        ToA_Matrix=(reshape(toa(1:tot_mesure),3,[]));
        
        ToA_Comp = [];
        
        ToA_Comp(1,:)=ToA_Matrix(1,:) - compensation_1(iii);
        ToA_Comp(2,:)=ToA_Matrix(2,:) - compensation_2(iii);
        ToA_Comp(3,:)=ToA_Matrix(3,:) - compensation_3(iii);
        
        toa1true=(sqrt(sum(([XR1,YR1]-pos_UE).^2,2)))./(inv(92.16e6*16).*3e8);
        toa2true=(sqrt(sum(([XR2,YR2]-pos_UE).^2,2)))./(inv(92.16e6*16).*3e8);
        toa3true=(sqrt(sum(([XR3,YR3]-pos_UE).^2,2)))./(inv(92.16e6*16).*3e8);
        toa1_error = [toa1_error (ToA_Comp(1,:) - toa1true)];
        toa2_error = [toa2_error (ToA_Comp(2,:) - toa2true)];
        toa3_error = [toa3_error (ToA_Comp(3,:) - toa3true)];
        
        DOA12=(ToA_Comp(1,:)-ToA_Comp(2,:) + 21.64).*inv(92.16e6*16).*3e8;
        DOA13=(ToA_Comp(1,:)-ToA_Comp(3,:) + 14.84).*inv(92.16e6*16).*3e8;
        DOA23=(ToA_Comp(2,:)-ToA_Comp(3,:) + 6.79).*inv(92.16e6*16).*3e8;
        
        doa12true=(sqrt(sum(([XR1,YR1]-pos_UE).^2,2))-sqrt(sum(([XR2,YR2]-pos_UE).^2,2)));
        doa13true=(sqrt(sum(([XR1,YR1]-pos_UE).^2,2))-sqrt(sum(([XR3,YR3]-pos_UE).^2,2)));
        doa23true=(sqrt(sum(([XR2,YR2]-pos_UE).^2,2))-sqrt(sum(([XR3,YR3]-pos_UE).^2,2)));
        
        doa12_error = [doa12_error (DOA12 - doa12true)];
        doa13_error = [doa13_error (DOA13 - doa13true)];
        doa23_error = [doa23_error (DOA23 - doa23true)];
        
        TDOA12=(ToA_Comp(1,:)-ToA_Comp(2,:));
        TDOA13=(ToA_Comp(1,:)-ToA_Comp(3,:));
        TDOA23=(ToA_Comp(2,:)-ToA_Comp(3,:));
        tdoa12true=(sqrt(sum(([XR1,YR1]-pos_UE).^2,2))-sqrt(sum(([XR2,YR2]-pos_UE).^2,2)))./(inv(92.16e6*16).*3e8);
        tdoa13true=(sqrt(sum(([XR1,YR1]-pos_UE).^2,2))-sqrt(sum(([XR3,YR3]-pos_UE).^2,2)))./(inv(92.16e6*16).*3e8);
        tdoa23true=(sqrt(sum(([XR2,YR2]-pos_UE).^2,2))-sqrt(sum(([XR3,YR3]-pos_UE).^2,2)))./(inv(92.16e6*16).*3e8);
        tdoa12_error = [tdoa12_error (TDOA12 - tdoa12true)];
        tdoa13_error = [tdoa13_error (TDOA13 - tdoa13true)];
        tdoa23_error = [tdoa23_error (TDOA23 - tdoa23true)];
        
        fprintf("POS: %d median(tdoa12_error) %.2f median(tdoa13_error) %.2f median(tdoa23_error) %.2f\n", ...
            reference_pos, median(TDOA12 - tdoa12true), median(TDOA13 - tdoa13true), median(TDOA23 - tdoa23true))
        
        F=find((~isnan(DOA12)).*(~isnan(DOA13)).*(~isnan(DOA23)));
        
        DOA12_w=movmean(DOA12(F),sliding_windows);
        DOA13_w=movmean(DOA13(F),sliding_windows);
        DOA23_w=movmean(DOA23(F),sliding_windows);
        
        [GX,GY]=meshgrid(0:0.2:12,0:0.2:6);
        GPz=[GX(:),GY(:)];
        
        pos_est = NaN(length(DOA12_w),2);
        
        for index_measure=1:length(DOA12_w)
            res = NaN(1,length(GPz));
            for index_pos_analize=1:length(GPz)
                d1=sqrt((XR1-GPz(index_pos_analize,1)).^2+(YR1-GPz(index_pos_analize,2)).^2);
                d2=sqrt((XR2-GPz(index_pos_analize,1)).^2+(YR2-GPz(index_pos_analize,2)).^2);
                d3=sqrt((XR3-GPz(index_pos_analize,1)).^2+(YR3-GPz(index_pos_analize,2)).^2);
                res(index_pos_analize)=((d1-d2)-DOA12_w(1,index_measure)).^2+...
                    ((d1-d3)-DOA13_w(1,index_measure)).^2+...
                    ((d2-d3)-DOA23_w(1,index_measure)).^2;
            end
            [~,i]=min(sum(res,1));
            pos_est(index_measure,:)=GPz(i,:);
        end
        
        for index_pos_analize=1:length(GPz)
            d1=sqrt((XR1-GPz(index_pos_analize,1)).^2+(YR1-GPz(index_pos_analize,2)).^2);
            d2=sqrt((XR2-GPz(index_pos_analize,1)).^2+(YR2-GPz(index_pos_analize,2)).^2);
            d3=sqrt((XR3-GPz(index_pos_analize,1)).^2+(YR3-GPz(index_pos_analize,2)).^2);
            res(index_pos_analize)=((d1-d2)-doa12true).^2+...
                ((d1-d3)-doa13true).^2+...
                ((d2-d3)-doa23true).^2;
        end
        [~,i]=min(sum(res,1));
        pos_est_ideal=GPz(i,:);
        
        pos_err_tmp=sqrt(sum((pos_UE-pos_est).^2,2)).';
        
        pos_err = [pos_err pos_err_tmp];
        
        if iii == 2 && reference_pos ~= 1 && reference_pos ~= 2 && reference_pos ~= 16
            figure(105)
            MEAN_ERROR = mean(pos_err_tmp);
            if reference_pos == 3
                plot(GPz(:,1),GPz(:,2),'*','linewidth',0.5);
                hold on
                plot(XR1,YR1,'^y','MarkerSize',8,'MarkerFaceColor','y')
                plot(XR2,YR2,'^y','MarkerSize',8,'MarkerFaceColor','y')
                plot(XR3,YR3,'^y','MarkerSize',8,'MarkerFaceColor','y')
                xlabel("x [m]")
                ylabel("y [m]")
                set(gca,'fontsize',24);
                xlim([-0.1 12.1])
                ylim([-0.1 6.1])
            end
            half_val_scala = 2;
            if MEAN_ERROR < half_val_scala
                red_value = MEAN_ERROR/half_val_scala;
                green_value = 1;
            else
                red_value = 1;
                green_value = max(0, 1 - ((MEAN_ERROR-half_val_scala)/half_val_scala));
            end
            plot(pos_UE(1),pos_UE(2),"square" , 'MarkerSize',15, 'MarkerEdgeColor',[red_value,green_value,0],'MarkerFaceColor',[red_value,green_value,0])
            
            if reference_pos == 17
                A = imread('gyrcb.png');
%                 imshow(A)
                % get base color table from image
                basect = im2double(permute(A(40,24:517,:),[2 3 1]));
                N = 256;
                nb = size(basect,1);
                % interpolate to get a specified-length version
                CT1 = interp1(1:nb,basect,linspace(1,nb,N));
                % also make a symmetric version of the same length
                CT2 = interp1(1:2*nb,[flipud(basect); basect],linspace(1,2*nb,N));
                % Apply the new CT
%                 surf(peaks)
                colormap(CT1*2)
                colorbar
            end
            
        end
        
%         if iii == 2
%             if reference_pos == 15 || reference_pos == 14
%                 figure;
%                 plot(GPz(:,1),GPz(:,2),'*','linewidth',0.5);
%                 hold on
%                 plot(pos_UE(1),pos_UE(2),'*g','linewidth',3);
%                 plot(movmean(pos_est(:,1),100),movmean(pos_est(:,2),100),'*r','linewidth',3);
%                 plot(XR1,YR1,'^y','MarkerSize',8,'MarkerFaceColor','y')
%                 plot(XR2,YR2,'^y','MarkerSize',8,'MarkerFaceColor','y')
%                 plot(XR3,YR3,'^y','MarkerSize',8,'MarkerFaceColor','y')
%                 xlabel("x [m]")
%                 ylabel("y [m]")
%                 set(gca,'fontsize',24);
%                 xlim([-0.1 12.1])
%                 ylim([-0.1 6.1])
%             end
%             
%         end
        
    end
    
    %---------- code for ECDF ----------%
    figure(101)
    PLOT = pos_err;
    PLOT(isnan(PLOT)) = [];
    if numel(PLOT) > 0
        [a,b] = ecdf(PLOT);
        plot(b,a,'linewidth',3);
        hold on
    end
        
end
    
grid on;
set(gca,'fontsize',24);
xlabel('Distance error [m]');
ylabel('ECDF');
legend("X","Y")
