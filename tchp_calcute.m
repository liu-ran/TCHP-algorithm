function [TCHP,d_t_ref]=tchp_calcute(t,s,lev,lat,t_ref);
%% t ,s (k,i,j)

[km,im,jm]=size(t);
d_t26=zeros(im,jm);
%º∆À„D26/28
for i=1:im
    for j=1:jm
        if isnan(t(1,i,j))==0
            if max(squeeze(t(:,i,j)))<t_ref
                d_t26(i,j)=nan;
            else
                if t(1,i,j)==t_ref
                    d_t26(i,j)=0;
                else
                    it=max(find(t(:,i,j)>t_ref));
                    if isnan(t(it+1,i,j))==1
                        d_t26(i,j)=lev(it);
                    else
                        T=it:it+1;
                        d_t26(i,j)=spline(t(T,i,j),lev(T),t_ref);
                    end
                end
            end
        else
            d_t26(i,j)=nan;
        end
    end
end
d_t_ref=d_t26;
%º∆À„TCHP
TCHP=zeros(im,jm);
d=0;
pres=0;
den=0;
Cp=0;
for i=1:im
    for j=1:jm
        if isnan(d_t26(i,j))==0
            d=d_t26(i,j);
            it=max(find(lev<d));
            for k=1:it
                pres=sw_pres(lev(k),lat(i));
                den=sw_dens(s(k,i,j),t(k,i,j),pres);
                %Cp=4180;
                Cp=sw_cp(s(k,i,j),t(k,i,j),pres);
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
                TCHP(i,j)=TCHP(i,j)+den*Cp*(t(k,i,j)-t_ref)*dz;
            end
        else
            TCHP(i,j)=nan;
        end
    end
end
