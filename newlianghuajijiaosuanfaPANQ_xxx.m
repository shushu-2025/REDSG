clc;        % clear command window
clear all;  % clear all variables
%This program tests the performance of the REDSG2 algorithm, 
%corresponding to Figure 6(a) in the article.
N = 15;             % number of agents in the network
L = 6;%64*2;             % size of wo, which is Mx1

wo = randn(L,N);wo_m= randn(L,N);%
for k=1:N
    wo_m(:,k) =1*wo(:,k)/norm(wo(:,k));  % real data (vector)
end
load wo_buton wo_m;
%
% sigma_out = rand(N,1)*30;
% save sigma_out sigma_out;
% RU= rand(1,N);
% save RU RU;
% load topology8;
% A=zeros(N,N);
% deg=sum(adjacency);
% for k=1:N
%     for l=1:N
%         if adjacency(l,k)~=0
%             if l==k
%                A(l,k)=1-(deg(k)-1)/N;
%             else
%                A(l,k)=1/N;
%             end
%         else
%             A(l,k)=0;
%         end
%     end
% end
% A_ba=(eye(N)+A)./2;
% save A_ba A_ba;

%load wo;
load sigma_rrandn%sigma_r8_12;%sigma_rnewt;%sigma_rnew2;
snrlow=40;snrup=60;%1e5*1*sigma_r;%
%snrlow=0;snrup=20;
sigma_r=(snrup-snrlow)*rand(N,N)+snrlow;
for ell=1:N
    sigma_r(ell,ell)=0;
end
sigma_rlh=1+rand(1,N);
load sigma_out;
load RU;
load topology15;
load A_ba;
nei_num=sum(adjacency,2);
p=1/N*ones(N,1);%nei_num/sum(nei_num);
sum1=0;
sum2=0;c=10;
for k=1:N
    sum1=sum1+p(k,1)*RU(1,k);
    sum2=sum2+p(k,1)*RU(1,k)*c*wo_m(:,k);
    %         sum1=sum1+mu_dlms*p(k,1)*R(1,k);
    %         sum2=sum2+mu_dlms*p(k,1)*R(1,k)*wo_m(:,k);
end
wo=inv(diag(sum1))*sum2;


mu=0.01;0.02;
%sigma_r=1e-10*ones(N ,N );
%  sigma_r=abs(1e-7*randn(N ,N ));%1e-12+(1e-8-1e-12)*rand(N ,N );
%  for i=1:N
%      sigma_r(i,i)=0;
%  end
%sigma_r(4,1)=0.0000005;
Num_trial = 50;
Num_iter  =100000;%10000000;%
MSD_agent= zeros(N,Num_iter);MSD_agent_c= zeros(N,Num_iter);
MSD_agent_dd= zeros(N,Num_iter);MSD_agent_ex= zeros(N,Num_iter);
MSD_agent_exTiC= zeros(N,Num_iter);MSD_agent_ced= zeros(N,Num_iter);
sav_RCDT2= zeros(1,Num_iter);sav_TiC=zeros(1,Num_iter);
sav_RCDT3=zeros(1,Num_iter);sav_ic=zeros(1,Num_iter);sav_P2=zeros(1,Num_iter);
sav_CED=zeros(1,Num_iter);
MSD_agent_cwQ2= zeros(N,Num_iter);MSD_agent_cwQ2_IC= zeros(N,Num_iter);
omiga=0.02;eta=1/(2*L);bit=2;
q_type = [12,6];
qc = quantizer('fixed','round','saturate',q_type);
for tt=1:Num_trial % iterating over experiments
    tt

    %     for k=1:N
    %         xk(:,k)=sqrt(RU(1,k)).*randn(Num_iter+L,1);
    %     end
    %     %%
    %     dk=zeros(Num_iter+L,N);
    %     for k=1:N
    %         y =  filter(1*wo_m(:,k),1,xk(:,k));%
    %         var_y=sum(y.^2)/(Num_iter+L);
    %         sigma_vo=1/(10^(sigma_out(k)/10));%var_y
    %         dk(:,k)=y+randn(Num_iter+L,1)*sqrt(100*sigma_vo);
    %     end

    %%
    cphi_y=zeros(L,N);cw=0*randn(L,N);
    w_dd=cw;y_dd=zeros(L,N);pphi_dd=zeros(L,N);%zeros(L,N)
    phi_pex=0*randn(L,N);wex=cw;yex=zeros(L,N);
    yexTiC=zeros(L,N);
    wexced=cw;yexced=zeros(L,N);
    wexTiC=cw;yexjq=zeros(L,N);
    cphi_yjq=zeros(L,N);cwjq=cw;
    sphiex=0*randn(L,N);uex=0*randn(L,N);
    w=cw;phi_yi=zeros(L,N);%w_p=zeros(L,N,(Num_iter+L));

    %%
    sum_linkav=0*randn(L,N);
    sum_link=0*randn(L,N);
    yy_dd=0*randn(L,N);
    wexp=0*randn(L,N);wexpp=0*randn(L,N);
    %%
    cphi_yQ2jq=0*randn(L,N);cphi_yQ2=0*randn(L,N);
    cphi_yQ2_IC=0*randn(L,N);
    prdettdu=zeros(L,N);
    prdettdu_IC=zeros(L,N);
    cwQ2_IC=0*randn(L,N);cwQ2=0*randn(L,N);
    isav_RCDT2= 0;isav_TiC=0;isav_CED=0;
    isav_RCDT3=0;isav_ic=0;isav_P2=0;
    %     w_la=w;
    %     w_lala=w_la;
    for i=1:(Num_iter)
        %%
        for k=1:N
            xk(:,k)=sqrt(RU(1,k)).*randn(L,1);
        end
        %%
        %dk=zeros(Num_iter+L,N);
        if i<=2.5e5%;3*Num_iter/2%
            for k=1:N
                y = xk(:,k)'*c*wo_m(:,k);%wo;% filter(,1,xk(:,k));%
                %var_y=sum(y.^2)/(Num_iter+L);
                sigma_vo=1/(10^(sigma_out(k)/5/10));
                dk(k)=y+randn(1,1)*sqrt(10*sigma_vo);
            end
        else
            for k=1:N
                y = xk(:,k)'*(-c)*wo_m(:,k);%wo;% filter(,1,xk(:,k));%
                %var_y=sum(y.^2)/(Num_iter+L);
                sigma_vo=1/(10^(sigma_out(k)/5/10));
                dk(k)=y+randn(1,1)*sqrt(10*sigma_vo);
            end
        end
        eps11=0.05;%0.005;%
        eps1=2*eps11;%0.05;%0.2*mu^2;%20/(1*i+200);%1/(0.1*i^0.6+1);
        eps2=0.9;%0.9;%0.8;%0.005
        %%
        if i<=2.5e5%3*Num_iter/2
            for k=1:N
                MSD_agent_c(k,i)=MSD_agent_c(k,i)+(wo-cw(:,k))'*(wo-cw(:,k));
                MSD_agent_exTiC(k,i)=MSD_agent_exTiC(k,i)+(wo-wexTiC(:,k))'*(wo-wexTiC(:,k));%TiC
                MSD_agent_ex(k,i)=MSD_agent_ex(k,i)+(wo-wex(:,k))'*(wo-wex(:,k));%propsed
                MSD_agent_cwQ2(k,i)=MSD_agent_cwQ2(k,i)+(wo-cwQ2(:,k))'*(wo-cwQ2(:,k));
                MSD_agent_cwQ2_IC(k,i)=MSD_agent_cwQ2_IC(k,i)+(wo-cwQ2_IC(:,k))'*(wo-cwQ2_IC(:,k));


                MSD_agent_ced(k,i)=MSD_agent_ced(k,i)+(wo-wexced(:,k))'*(wo-wexced(:,k));%CED
            end
        else
            for k=1:N
                MSD_agent_c(k,i)=MSD_agent_c(k,i)+(-wo-cw(:,k))'*(-wo-cw(:,k));
                MSD_agent_exTiC(k,i)=MSD_agent_exTiC(k,i)+(-wo-wexTiC(:,k))'*(-wo-wexTiC(:,k));%TiC
                MSD_agent_ex(k,i)=MSD_agent_ex(k,i)+(-wo-wex(:,k))'*(-wo-wex(:,k));%propsed
                MSD_agent_cwQ2(k,i)=MSD_agent_cwQ2(k,i)+(-wo-cwQ2(:,k))'*(-wo-cwQ2(:,k));
                MSD_agent_cwQ2_IC(k,i)=MSD_agent_cwQ2_IC(k,i)+(-wo-cwQ2_IC(:,k))'*(-wo-cwQ2_IC(:,k));


                MSD_agent_ced(k,i)=MSD_agent_ced(k,i)+(-wo-wexced(:,k))'*(-wo-wexced(:,k));%CED

            end
        end
        %         w_lala=w_la;
        %         w_la=w;
        %         for k=1:N
        %             phi_w(:,k)=w(:,k)+mu*xk(i:-1:i-L+1,k)*(dk(i,k)-xk(i:-1:i-L+1,k)'*w(:,k));
        %         end
        %         phi_yp=zeros(L,N);
        %         for k=1:N
        %             for ell=1:N
        %                 if A_ba(ell,k)~=0
        %                     phi_yp(:,k)=phi_yp(:,k)+A_ba(ell,k)*(phi_w(:,ell)+sqrt(sigma_r(ell,k))*randn(L,1));
        %                 end
        %             end
        %         end
        %         w=phi_yp;
        %%         %%  %%%%%%%%%%%%%%%%%%%%%% c-GT
        if mod(i,2)==0
            %%  %%%%%%%%%%%%%%%%%%%%%% c-SGT
            for k=1:N
                [rnoiseex2,rx1(k)]=PANQ(cphi_y(:,k),omiga,eta,L);
                noiseex2=rnoiseex2-cphi_y(:,k);
                %noiseex=rq(wex(:,k),3*bit)-wex(:,k); %
                %noiseex2 =quantize(qc,shphiex(:,k))-shphiex(:,k);
                %noiseex=sqrt(1e-3*sigma_rlh(1,k))*randn(L,1);%
                cphi_ylihua(:,k)=cphi_y(:,k)+noiseex2;%norm(psiex(:,k)-1*wo)^2
            end
            cphi_yp=zeros(L,N);
            for k=1:N
                for ell=1:N
                    if A_ba(ell,k)~=0
                        cphi_yp(:,k)=cphi_yp(:,k)+A_ba(ell,k)*(cphi_ylihua(:,ell));
                    end
                end
            end
            for k=1:N
                cphi_yp(:,k)=cphi_yp(:,k)-cphi_ylihua(:,k)+cphi_y(:,k);
            end
            %
            for k=1:N
                cphi_yz(:,k)=(1-eps2*1)*cphi_y(:,k)+eps2*1*cphi_yp(:,k)-(1-0)*xk(:,k)*(dk(k)-xk(:,k)'*cw(:,k));%*(RU(1,k)*eye(L))*(10*wo_m(:,k)-cw(:,k));%
                %;%
                %(RU(1,k)*eye(L))*(1*wo_m(:,k)-cw(:,k));
            end
            for k=1:N
                [rnoiseex2,rx2(k)]=PANQ(cw(:,k),omiga,eta,L);
                noiseex2 =rnoiseex2-cw(:,k);
                %noiseex=rq(wex(:,k),3*bit)-wex(:,k); %
                %noiseex2 =quantize(qc,shphiex(:,k))-shphiex(:,k);
                %noiseex=sqrt(1e-3*sigma_rlh(1,k))*randn(L,1);%
                cwlihua(:,k)=cw(:,k)+noiseex2;%norm(psiex(:,k)-1*wo)^2
            end
            isav_RCDT2=isav_RCDT2+sum(rx1+rx2);
            sav_RCDT2(i)=isav_RCDT2;
            cphi_wp=zeros(L,N);
            for k=1:N
                for ell=1:N
                    if A_ba(ell,k)~=0%
                        cphi_wp(:,k)=cphi_wp(:,k)+A_ba(ell,k)*(cwlihua(:,ell));
                    end
                end
            end
            for k=1:N
                cphi_wp(:,k)=cphi_wp(:,k)-cwlihua(:,k)+cw(:,k);
            end
            for k=1:N
                %cw(:,k)= cw(:,k)-mu*eps2*(cphi_yz(:,k)-1*cphi_y(:,k));%SGT
                cw(:,k)= (1-eps1)*cw(:,k)+eps1*cphi_wp(:,k)-mu*eps1*(cphi_yz(:,k)-1*cphi_y(:,k));%GT
                %cw(:,k)= cphi_wp(:,k)-mu*eps2*(cphi_yz(:,k)-1*cphi_y(:,k));
            end
            cphi_y=cphi_yz;
            %% RCSGT3
            for k=1:N
                [rnoiseex2,rx1(k)]=PANQ(cphi_yQ2(:,k),omiga,eta,L);
                noiseex2=rnoiseex2-cphi_yQ2(:,k);
                cphi_ylihuaQ2(:,k)=cphi_yQ2(:,k)+noiseex2;%norm(psiex(:,k)-1*wo)^2
            end
            cphi_ypQ2=zeros(L,N);
            for k=1:N
                for ell=1:N
                    if A_ba(ell,k)~=0
                        cphi_ypQ2(:,k)=cphi_ypQ2(:,k)+A_ba(ell,k)*(cphi_ylihuaQ2(:,ell));
                    end
                end
            end
            for k=1:N
                cphi_ypQ2(:,k)=cphi_ypQ2(:,k)-cphi_ylihuaQ2(:,k)+cphi_yQ2(:,k);
            end
            for k=1:N
                dettdu(:,k)=-(1-0)*xk(:,k)*(dk(k)-xk(:,k)'*cwQ2(:,k));
                cphi_yQ2(:,k)=(1-eps2*1)*cphi_yQ2(:,k)+eps2*1*cphi_ypQ2(:,k)+dettdu(:,k)-prdettdu(:,k);%*(RU(1,k)*eye(L))*(10*wo_m(:,k)-cw(:,k));%
            end
            prdettdu=dettdu;
            for k=1:N
                [rnoiseex2,rx2(k)]=PANQ(cwQ2(:,k),omiga,eta,L);
                noiseex2=rnoiseex2-cwQ2(:,k);
                cwlihuaQ2(:,k)=cwQ2(:,k)+noiseex2;%norm(psiex(:,k)-1*wo)^2
            end
            isav_RCDT3=isav_RCDT3+sum(rx1+rx2);
            sav_RCDT3(i)=isav_RCDT3;
            cphi_wpQ2=zeros(L,N);
            for k=1:N
                for ell=1:N
                    if A_ba(ell,k)~=0%
                        cphi_wpQ2(:,k)=cphi_wpQ2(:,k)+A_ba(ell,k)*(cwlihuaQ2(:,ell));
                    end
                end
            end
            for k=1:N
                cphi_wpQ2(:,k)=cphi_wpQ2(:,k)-cwlihuaQ2(:,k)+cwQ2(:,k);
            end
            for k=1:N
                cwQ2(:,k)= (1-eps1)*cwQ2(:,k)+eps1*cphi_wpQ2(:,k)-mu*eps1*cphi_yQ2(:,k);%GT
            end
            %%
            %% _IC
            for k=1:N
                [rnoiseex2,rx1(k)]=PANQ(cphi_yQ2_IC(:,k),omiga,eta,L);
                noiseex2=rnoiseex2-cphi_yQ2_IC(:,k);
                cphi_ylihuaQ2_IC(:,k)=cphi_yQ2_IC(:,k)+noiseex2;%norm(psiex(:,k)-1*wo)^2
            end
            cphi_ypQ2_IC=zeros(L,N);
            for k=1:N
                for ell=1:N
                    if A_ba(ell,k)~=0
                        cphi_ypQ2_IC(:,k)=cphi_ypQ2_IC(:,k)+A_ba(ell,k)*(cphi_ylihuaQ2_IC(:,ell));
                    end
                end
            end
            for k=1:N
                dettdu_IC(:,k)=-(1-0)*xk(:,k)*(dk(k)-xk(:,k)'*cwQ2_IC(:,k));
                cphi_yQ2_IC(:,k)=(1-eps1*1)*cphi_yQ2_IC(:,k)+eps1*1*cphi_ypQ2_IC(:,k)+dettdu_IC(:,k)-prdettdu_IC(:,k);%*(RU(1,k)*eye(L))*(10*wo_m(:,k)-cw(:,k));%
            end
            prdettdu_IC=dettdu_IC;
            for k=1:N
                [rnoiseex2,rx2(k)]=PANQ(cwQ2_IC(:,k),omiga,eta,L);
                noiseex2=rnoiseex2-cwQ2_IC(:,k);
                cwlihuaQ2_IC(:,k)=cwQ2_IC(:,k)+noiseex2;%norm(psiex(:,k)-1*wo)^2
            end
            isav_ic=isav_ic+sum(rx1+rx2);
            sav_ic(i)=isav_ic;
            cphi_wpQ2_IC=zeros(L,N);
            for k=1:N
                for ell=1:N
                    if A_ba(ell,k)~=0%
                        cphi_wpQ2_IC(:,k)=cphi_wpQ2_IC(:,k)+A_ba(ell,k)*(cwlihuaQ2_IC(:,ell));
                    end
                end
            end
            for k=1:N
                cwQ2_IC(:,k)= (1-eps1)*cwQ2_IC(:,k)+eps1*cphi_wpQ2_IC(:,k)-mu*eps1*cphi_yQ2_IC(:,k);%GT
            end
        end
        %%
        %% PP 共享一个变量
        % %         for k=1:N
        % %             phi_dd(:,k)=w_dd(:,k)+mu*xk(i:-1:i-L+1,k)*(dk(i,k)-xk(i:-1:i-L+1,k)'*w_dd(:,k))-y_dd(:,k);
        % %         end
        % %         cphi_dd=zeros(L,N);
        % %         for k=1:N
        % %             for ell=1:N
        % %                 if A_ba(ell,k)~=0
        % %                     cphi_dd(:,k)=cphi_dd(:,k)+A_ba(ell,k)*(phi_dd(:,ell)+sqrt(1*sigma_r(ell,k))*randn(L,1));
        % %                 end
        % %             end
        % %         end
        % %         for k=1:N
        % %             w_dd(:,k)=(1-eps1)*w_dd(:,k)+eps1*cphi_dd(:,k);
        % %         end
        % %         for k=1:N
        % %             y_dd(:,k)=y_dd(:,k)+eps1*(phi_dd(:,k)-cphi_dd(:,k));%phi_dd(:,k)-w_dd(:,k);%
        % %         end
        %%
        %%
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Propo
        %A_baex=(1-(1))*eye(N)+(1)*A_ba;%0.0001
        for k=1:N
            [rnoiseex2,rx1(k)]=PANQ(wex(:,k),omiga,eta,L);
            noiseex2=rnoiseex2-wex(:,k);
            %noiseex=rq(wex(:,k),3*bit)-wex(:,k); %
            %noiseex2 =quantize(qc,shphiex(:,k))-shphiex(:,k);
            %noiseex=sqrt(1e-3*sigma_rlh(1,k))*randn(L,1);%
            wexlihua(:,k)=wex(:,k)+noiseex2;%norm(psiex(:,k)-1*wo)^2
        end
        isav_P2=isav_P2+sum(rx1);
        sav_P2(i)=isav_P2;
        cphi_ex=zeros(L,N);
        for k=1:N
            for ell=1:N
                if A_ba(ell,k)~=0
                    cphi_ex(:,k)=cphi_ex(:,k)+A_ba(ell,k)*(wexlihua(:,ell));
                end
            end
        end
        for k=1:N
            cphi_ex(:,k)=cphi_ex(:,k)-wexlihua(:,k)+wex(:,k);
        end
        zhjg=cphi_ex;%pyex=yex;
        for k=1:N
            yex(:,k)=yex(:,k)+eps11*(wex(:,k)-zhjg(:,k));
            wex(:,k)=(1-eps11)*wex(:,k)+eps11*zhjg(:,k)+eps11*mu*xk(:,k)*(dk(k)-xk(:,k)'*wex(:,k))...*(RU(1,k)*eye(L))*(10*wo_m(:,k)-wex(:,k))...%
                -eps11*yex(:,k);
        end
        %% TiC
        for k=1:N
            [rnoiseex2,rx1(k)]=PANQ(wexTiC(:,k),omiga,eta,L);
            noiseex2=rnoiseex2-wexTiC(:,k);
            %noiseex=rq(wex(:,k),3*bit)-wex(:,k); %
            %noiseex2 =quantize(qc,shphiex(:,k))-shphiex(:,k);
            %noiseex=sqrt(1e-3*sigma_rlh(1,k))*randn(L,1);%
            wexlihuaTiC(:,k)=wexTiC(:,k)+noiseex2;%norm(psiex(:,k)-1*wo)^2
        end
        isav_TiC=isav_TiC+sum(rx1);
        sav_TiC(i)=isav_TiC;
        cphi_ex=zeros(L,N);
        for k=1:N
            for ell=1:N
                if A_ba(ell,k)~=0
                    cphi_ex(:,k)=cphi_ex(:,k)+(wexlihuaTiC(:,ell));
                end
            end
        end

        
        for k=1:N
            Rsum(:,k)=cphi_ex(:,k);
            zhjgTiC(:,k)=sum(sign(A_ba(:,k)))*wexlihuaTiC(:,k)-Rsum(:,k);
        end
        for k=1:N
            wexTiC(:,k)=0.99*wexTiC(:,k)+(1-0.99)*wexlihuaTiC(:,k)-eps11*mu*(100*zhjgTiC(:,k)-xk(:,k)*(dk(k)-xk(:,k)'*wexTiC(:,k))...*(RU(1,k)*eye(L))*(10*wo_m(:,k)-wex(:,k))...%
                +yexTiC(:,k));
            yexTiC(:,k)=yexTiC(:,k)+eps11*zhjgTiC(:,k);%eps11
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CEDAS
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %A_baex=(1-(1))*eye(N)+(1)*A_ba;%0.0001
        for k=1:N
            cwexced(:,k)=wexced(:,k)+eps11*mu*xk(:,k)*(dk(k)-xk(:,k)'*wexced(:,k))...*(RU(1,k)*eye(L))*(10*wo_m(:,k)-wex(:,k))...%
                -yexced(:,k);
        end

        for k=1:N
            [rnoiseex2,rx1(k)]=PANQ(cwexced(:,k),omiga,eta,L);
            noiseex2=rnoiseex2-cwexced(:,k);
            %noiseex=rq(wex(:,k),3*bit)-wex(:,k); %
            %noiseex2 =quantize(qc,shphiex(:,k))-shphiex(:,k);
            %noiseex=sqrt(1e-3*sigma_rlh(1,k))*randn(L,1);%
            cwexlihuaced(:,k)=cwexced(:,k)+noiseex2;%norm(psiex(:,k)-1*wo)^2
        end
        isav_CED=isav_CED+sum(rx1);
        sav_CED(i)=isav_CED;
        cphi_ex=zeros(L,N);
        for k=1:N
            for ell=1:N
                if A_ba(ell,k)~=0
                    cphi_ex(:,k)=cphi_ex(:,k)+A_ba(ell,k)*(cwexlihuaced(:,ell));
                end
            end
        end
        for k=1:N
            cphi_ex(:,k)=cphi_ex(:,k)-cwexlihuaced(:,k)+cwexced(:,k);
        end
        czhjg=cphi_ex;%pyex=yex;
        for k=1:N
            yexced(:,k)=yexced(:,k)+eps11/2*(cwexced(:,k)-czhjg(:,k));
            wexced(:,k)=wexced(:,k)+eps11*mu*xk(:,k)*(dk(k)-xk(:,k)'*wexced(:,k))...*(RU(1,k)*eye(L))*(10*wo_m(:,k)-wex(:,k))...%
                -yexced(:,k);
        end
    end
end

tt2c_ced = MSD_agent_ced/Num_trial; % each row contains the MSD evolution of the corresponding agent
MSD_av2c_ced = sum(tt2c_ced)/N;   % learning curve for adaptive combination matrix averaged over experiments and number of nodes
MSD_av_dbc_ced = 10*log10(MSD_av2c_ced);

%%
tt2c = MSD_agent_c/Num_trial; % each row contains the MSD evolution of the corresponding agent
MSD_av2c = sum(tt2c)/N;   % learning curve for adaptive combination matrix averaged over experiments and number of nodes
MSD_av_dbc = 10*log10(MSD_av2c);
%%

tt2ex = MSD_agent_ex/Num_trial; % each row contains the MSD evolution of the corresponding agent
MSD_av2ex = sum(tt2ex)/N;   % learning curve for adaptive combination matrix averaged over experiments and number of nodes
MSD_av_dbex = 10*log10(MSD_av2ex);%filter( 1/1000*ones(1000,1),1, MSD_av_dbex(1,L+1:end))

tt2exTiC = MSD_agent_exTiC/Num_trial; % each row contains the MSD evolution of the corresponding agent
MSD_av2exTiC = sum(tt2exTiC)/N;   % learning curve for adaptive combination matrix averaged over experiments and number of nodes
MSD_av_dbexTiC = 10*log10(MSD_av2exTiC);

tt2cwQ2 = MSD_agent_cwQ2/Num_trial; % each row contains the MSD evolution of the corresponding agent
MSD_av2cwQ2 = sum(tt2cwQ2)/N;   % learning curve for adaptive combination matrix averaged over experiments and number of nodes
MSD_av_dbcwQ2 = 10*log10(MSD_av2cwQ2);

tt2cwQ2_IC = MSD_agent_cwQ2_IC/Num_trial; % each row contains the MSD evolution of the corresponding agent
MSD_av2cwQ2_IC = sum(tt2cwQ2_IC)/N;   % learning curve for adaptive combination matrix averaged over experiments and number of nodes
MSD_av_dbcwQ2_IC = 10*log10(MSD_av2cwQ2_IC);

figure
hold on
iter = 1:Num_iter;
p1=plot(iter,MSD_av_dbc(1,1:end),'-.x');p1.MarkerSize = 10;
p1.MarkerIndices = 1:floor(Num_iter/3):length(iter);

p2=plot(iter,MSD_av_dbcwQ2(1,1:end),'-.x');p2.MarkerSize = 10;
p2.MarkerIndices = 1:floor(Num_iter/3):length(iter);

p3=plot(iter,MSD_av_dbcwQ2_IC(1,1:end),'-.x');p3.MarkerSize = 10;
p3.MarkerIndices = 1:floor(Num_iter/3):length(iter);

p4=plot(iter,MSD_av_dbc_ced(1,1:end),'-.x');p4.MarkerSize = 10;
p4.MarkerIndices = 1:floor(Num_iter/3):length(iter);

p5=plot(iter,MSD_av_dbexTiC(1,1:end),'-.x');p5.MarkerSize = 10;
p5.MarkerIndices = 1:floor(Num_iter/4):length(iter);

p6=plot(iter,MSD_av_dbex(1,1:end),'r-.x');p6.MarkerSize = 10;
p6.MarkerIndices = 1:floor(Num_iter/3):length(iter);
% % figure
% % iter = 1:Num_iter;
% % plot(iter, MSD_av_dbc(1,1:end),'b',...
% %     iter,MSD_av_dbcwQ2(1,1:end),'c',...
% %     iter,MSD_av_dbcwQ2_IC(1,1:end),'k',...
% %     iter,MSD_av_dbc_ced(1,1:end),'m',...
% %     iter,MSD_av_dbexTiC(1,1:end),'g',...
% %     iter,MSD_av_dbex(1,1:end),'r');
legend('RC-SGT2','RC-SGT3','IC-GT','CEDAS','TiCoPD','REDSG2');
xlabel('Rounds of communication');
ylabel('MSD(dB)');

figure
hold on
% % p111=plot(sav_RCDT2(2:2:end)/Num_trial,MSD_av_dbc(2:2:end),'k-.s');p111.MarkerSize = 10;
% % p111.MarkerIndices = 1:Num_iter/10:length(iter);%保护DLMS
% %
plot(sav_RCDT2(2:2:end)/Num_trial,MSD_av_dbc(2:2:end));plot(sav_RCDT3(2:2:end)/Num_trial,MSD_av_dbcwQ2(2:2:end));
plot(sav_ic(2:2:end)/Num_trial,MSD_av_dbcwQ2_IC(2:2:end));
plot(sav_CED(1:1:end)/Num_trial,MSD_av_dbc_ced(1:1:end));
plot(sav_TiC(1:1:end)/Num_trial,MSD_av_dbexTiC(1:1:end));
plot(sav_P2(1:1:end)/Num_trial,MSD_av_dbex(1,1:1:end));
legend('RC-SGT2','RC-SGT3','IC-GT','CEDAS','TiCoPD ','REDSG2');

%save data20240607






