clc;clear all;close all;
r=6378155;
nameall1=dir('G:\CMIP5_DATA\NEMO\T_allrun\');
nameall2=dir('G:\CMIP5_DATA\NEMO\S_allrun\');
lon=nc_varget(strcat('G:\CMIP5_DATA\NEMO\T_allrun\',nameall1(3).name),'lon');
lat=nc_varget(strcat('G:\CMIP5_DATA\NEMO\T_allrun\',nameall1(3).name),'lat');
lev=nc_varget(strcat('G:\CMIP5_DATA\NEMO\T_allrun\',nameall1(3).name),'deptht');
lat=lat(30:359);
nn1={'control run','WIND+SST+SSS','WIND','SST','SSS' };
nn2={'WIND+SST+SSS','WIND','SST','SSS'};
t_old_run0=nc_varget(strcat('G:\CMIP5_DATA\NEMO\T_allrun\',nameall1(3).name),'temp');
t_run0=t_old_run0(:,30:359,:);
d_t26=zeros(330,720);
for i=1:330
    for j=1:720
        if isnan(t_run0(1,i,j))==0
            if max(squeeze(t_run0(:,i,j)))<26
                d_t26(i,j)=nan;
            else
                if t_run0(1,i,j)==26
                    d_t26(i,j)=0;
                else
                    it=max(find(t_run0(:,i,j)>26));
                    if isnan(t_run0(it+1,i,j))==1
                        d_t26(i,j)=lev(it);
                    else
                        T=it:it+1;
                        d_t26(i,j)=spline(t_run0(T,i,j),lev(T),26);
                    end
                end
            end
        else
            d_t26(i,j)=nan;
        end
    end
end
figure(7)
set(gcf,'color','w')
[X3,Y3]=meshgrid(lon,lat);
m_proj('miller','long',[0.25 359.75],'lat',[-40.25 40.25]);
% m_gshhs_i('color','k','linewi',0.5);
C=m_contourf(X3,Y3,d_t26,[0:10:140]);
caxis([0,140]);
hold on
m_coast('patch',[0.72 0.72 0.72],'edgecolor','none');
m_grid('linewi',2,'linest','none','tickdir','out','fontsize',10);
colorbar;
set(colorbar,'Fontsize',14);
% hold on
% C= m_contour(X3,Y3,squeeze(ave2(28,71:232,241:541)),10,'w-','linewidth',1.5);
% clabel(C,'Fontsize',12);
title('D-T26-run0','Fontsize',14);
print(7,'-dpdf','D-T26-run0');
print(7,'-dtiff','D-T26-run0');
clear t_old_run0 t_run0

TCHP=zeros(330,720,5);
for n=3:7
    t_old=nc_varget(strcat('G:\CMIP5_DATA\NEMO\T_allrun\',nameall1(n).name),'temp');
    t=t_old(:,30:359,:);
    s_old=nc_varget(strcat('G:\CMIP5_DATA\NEMO\S_allrun\',nameall2(n).name),'salt');
    s=s_old(:,30:359,:);
    %计算TCHP_hisrcp
    
    d=0;
    pres=0;
    den=0;
    Cp=0;
    for i=1:330
        for j=1:720
            if isnan(d_t26(i,j))==0
                d=d_t26(i,j);
                it=max(find(lev<d));
                for k=1:it
                    pres=sw_pres(lev(k),lat(i));
                    den=sw_dens(s(k,i,j),t(k,i,j),pres);
                    Cp=4180;
                    %                 Cp=sw_cp(s(j,i,k),t(j,i,k),pres);
                    if it>2
                        if k==1
                            dz=(lev(2)-lev(1))/2;
                        elseif k==it
                            dz=(lev(k)-lev(k-1))/2+(d-lev(k));
                        elseif k>1 & k<it
                            dz=(lev(k+1)-lev(k-1))/2;
                        end
                    elseif it==2
                        if k==1
                            dz=(lev(2)-lev(1))/2;
                        elseif k==it
                            dz=(lev(k)-lev(k-1))/2+(d-lev(k));
                        end
                    elseif it==1
                        dz=d;
                    end
                    TCHP(i,j,n-2)=TCHP(i,j,n-2)+den*Cp*(t(k,i,j)-26)*dz;
                end
            else
                TCHP(i,j,n-2)=nan;
            end
        end
    end
end
for n=3:7
    figure(1)
    subplot(5,1,n-2)
    set(gcf,'color','w')
    [X3,Y3]=meshgrid(lon,lat);
    m_proj('miller','long',[0.25 359.75],'lat',[-40.25 40.25]);
    m_contourf(X3,Y3,squeeze(TCHP(:,:,n-2))/1e8,[0:1.5:12]);
    caxis([0,12]);
    m_coast('patch',[0.72 0.72 0.72],'edgecolor','none');
    m_grid('linewi',2,'linest','none','tickdir','out','fontsize',10);
    colorbar([0.75 0.1 0.02 0.82]);
    %     set(colorbar,'Fontsize',6);
    text(-2,0.6,nn1{n-2},'fontsize',8);
    %加最外围等值线
    if n>3
        tt=squeeze(TCHP(:,:,1));
        tt(isnan(tt))=0;
        hold on;
        C=m_contour(X3,Y3,tt,[0.0001 0.0001],'w-','linewidth',1.2);
        %         clabel(C,'manual','Fontsize',10)
    end
    if n==3
        title('NEMO TCHP ON THE SAME D26','Fontsize',14);
    end
    if n>3
        ttt2=squeeze(TCHP(:,:,n-2));
        ttt1=squeeze(TCHP(:,:,1));
        for i=1:330
            for j=1:720
                if isnan(ttt1(i,j))==0
                    dttt(i,j)=ttt2(i,j)-ttt1(i,j);
                else
                    dttt(i,j)=nan;
                end
            end
        end
        figure(2)
        subplot(4,1,n-3)
        set(gcf,'color','w')
        [X3,Y3]=meshgrid(lon,lat);
        m_proj('miller','long',[0.25 359.75],'lat',[-40.25 40.25]);
        m_contourf(X3,Y3,dttt/1e8,[-6:0.5:6]);
        caxis([-6,6]);
        m_coast('patch',[0.72 0.72 0.72],'edgecolor','none');
        m_grid('linewi',2,'linest','none','tickdir','out','fontsize',10);
        colorbar([0.8 0.1 0.02 0.82]);
        %         set(colorbar,'Fontsize',6);
        text(-2,0.6,nn2{n-3},'fontsize',8);
        if n>3
            hold on;
            C=m_contour(X3,Y3,dttt,[0.0001 0.0001],'r-','linewidth',1.2);
            %         clabel(C,'manual','Fontsize',10)
        end
        if n==4
            title('NEMO TCHP CHANGE ON THE SAME D26','Fontsize',14);
        end
    end
end
print(1,'-dpdf','NEMO_tchp_hisrcp');
print(1,'-dtiff','NEMO_tchp_hisrcp');
print(2,'-dpdf','NEMO_tchp_change_hisrcp');
print(2,'-dtiff','NEMO_tchp_change_hisrcp');


% %%%D26以上平均温度
% Ta=zeros(180,360);
% d=0;
% pres=0;
% den=0;
% Cp=0;
% for i=1:180
%     for j=1:360
%         if isnan(d_t26(i,j))==0
%             d=d_t26(i,j);
%             it=max(find(lev<d));
%             for k=1:it
%                 pres=sw_pres(lev(k),lat(i));
%                 den=sw_dens(s(j,i,k),t(j,i,k),pres);
%                 Cp=4180;
% %                 Cp=sw_cp(s(j,i,k),t(j,i,k),pres);
%                 if it>2
%                     if k==1
%                         dz=(lev(2)-lev(1))/2;
%                     elseif k==it
%                         dz=(lev(k)-lev(k-1))/2+(d-lev(k));
%                     elseif k>1 & k<it
%                         dz=(lev(k+1)-lev(k-1))/2;
%                     end
%                 elseif it==2
%                     if k==1
%                         dz=(lev(2)-lev(1))/2;
%                     elseif k==it
%                         dz=(lev(k)-lev(k-1))/2+(d-lev(k));
%                     end
%                 elseif it==1
%                     dz=d;
%                 end
%                 Ta(i,j)=Ta(i,j)+t(j,i,k)*dz/d_t26(i,j);
%             end
%         else
%             Ta(i,j)=nan;
%         end
%     end
% end
% figure(7)
% set(gcf,'color','w')
% [X3,Y3]=meshgrid(lon,lat);
% m_proj('miller','long',[0.5 359.5],'lat',[-40.5 40.5]);
% % m_gshhs_i('color','k','linewi',0.5);
% C=m_contourf(X3,Y3,Ta,[26:0.2:29]);
% % clabel(C,'manual','Fontsize',12)
% % m_pcolor(X3,Y3,dt3(:,:))
% % shading interp
% caxis([26,29]);
% hold on
% m_coast('patch',[0.72 0.72 0.72],'edgecolor','none');
% m_grid('linewi',2,'linest','none','tickdir','out','fontsize',10);
% colorbar;
% set(colorbar,'Fontsize',8);
% % hold on
% % C= m_contour(X3,Y3,squeeze(ave2(28,71:232,241:541)),10,'w-','linewidth',1.5);
% % clabel(C,'Fontsize',12);
% title('Ensemble mean temperature above D26','Fontsize',14);
% print('-dpdf','Ta-woa');
% print('-dtiff','Ta-woa');

% %面积平均
% r=6378155;
% tchp_MDR_woa=0;
% area_horizave=0;
% for i=100:116
%     for j=122:170
%         if isnan(TCHP(i,j))==0
%             area=(pi*r*cosd(lat(i))/180)*(pi*r/180);
%             area_horizave=area_horizave+area;
%             tchp_MDR_woa=tchp_MDR_woa+TCHP(i,j)*area;
%         end
%     end
% end
% tchp_MDR_woa=tchp_MDR_woa/area_horizave;
% 
% save TCHP.mat TCHP d_t26
% figure(2)
% set(gcf,'color','w')
% [X3,Y3]=meshgrid(lon,lat);
% m_proj('miller','long',[0.5 359.5],'lat',[-89.5 89.5]);
% % m_gshhs_i('color','k','linewi',0.5);
% C=m_contourf(X3,Y3,TCHP/1e8,25);
% clabel(C,'manual','Fontsize',12)
% % m_pcolor(X3,Y3,dt3(:,:))
% % shading interp
% caxis([0,10]);
% hold on
% m_coast('patch',[0.3 0.3 0.3],'edgecolor','none');
% m_grid('linewi',2,'linest','none','tickdir','out','fontsize',14);
% colorbar;
% set(colorbar,'Fontsize',18);
% % hold on
% % C= m_contour(X3,Y3,squeeze(ave2(28,71:232,241:541)),10,'w-','linewidth',1.5);
% % clabel(C,'Fontsize',12);
% title('TCHP-woa','Fontsize',18);
% print( '-djpeg', 'TCHP-woa');
% 
% %%%%%sum_TCHP
% r=6378155;
% sum_TCHP=0;
% for i=1:180
%     for j=1:360
%         if isnan(TCHP(i,j))==0
%             area=(pi*r*cosd(lat(i))/180)*(pi*r/180);
%             sum_TCHP=sum_TCHP+TCHP(i,j)*area;
%         end
%     end
% end
% 
% %%绘图
% figure(1)
% set(gcf,'color','w')
% [X3,Y3]=meshgrid(lon,lat);
% m_proj('miller','long',[0.5 359.5],'lat',[-39.5 39.5]);
% % m_gshhs_i('color','k','linewi',0.5);
% m_contourf(X3,Y3,squeeze(d_t26(:,:)),[0:10:140]);
% % clabel(C,'manual','Fontsize',12)
% caxis([0,140]);
% m_coast('patch',[0.72 0.72 0.72],'edgecolor','none');
% m_grid('linewi',2,'linest','none','tickdir','out','fontsize',10);
% colorbar;
% set(colorbar,'Fontsize',8);
% % hold on
% % C= m_contour(X3,Y3,squeeze(ave2(28,71:232,241:541)),10,'w-','linewidth',1.5);
% % clabel(C,'Fontsize',12);
% title('D26-woa','Fontsize',14);
% print( '-dpdf', 'D-T26-woa');
% print( '-dtiff', 'D-T26-woa');
% 
% figure(2)
% set(gcf,'color','w')
% [X3,Y3]=meshgrid(lon,lat);
% m_proj('miller','long',[0.5 359.5],'lat',[-39.5 39.5]);
% % m_gshhs_i('color','k','linewi',0.5);
% m_contourf(X3,Y3,TCHP/1e8,[0:1.5:12]);
% % clabel(C,'manual','Fontsize',12)
% caxis([0,12]);
% m_coast('patch',[0.72 0.72 0.72],'edgecolor','none');
% m_grid('linewi',2,'linest','none','tickdir','out','fontsize',10);
% colorbar;
% set(colorbar,'Fontsize',8);
% title('TCHP-woa','Fontsize',14);
% print( '-dpdf', 'TCHP-woa');
% print( '-dtiff', 'TCHP-woa');