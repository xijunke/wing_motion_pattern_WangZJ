%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wing_pattern_WangZJ
% ������������ֵ������Ӭ�������ŵĳ���˶�ѧ����
f=234;  % Hz����f��[0,inf];     % T=1/f;                            %��Frequency������ֵ            
phi_m=pi/2;                % rad.����phi_m��[0,pi/2];            %��Azimuthal amplitude������ֵ
% eta_m=72.7*pi/180;    % rad.����eta_m��[0,pi];                %��Pitching amplitude������ֵ
eta_m=45*pi/180;   
% eta_0=pi/2;             % rad.����eta_0��[eta_m-pi,pi-eta_m];     %��Pitching offset������ֵ
eta_0=0;                     % rad.������������������������������������Pitching offset��������ʼֵ����
% K=0.704;                     % K��[0,1]; % ��Affects the shape of phi(t )������ֵ
% C_eta=2.375;               % C_eta��[0,inf];                           %��Affects the duration of wing rotation������ֵ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phi_eta=-72.4*pi/180;  % rad.����Phi_eta��[-pi,pi];          %��Pitching phase offset������ֵ
% Phi_eta=0; 
% Phi_eta=pi;
% Phi_eta=pi/4; 
% Phi_eta=-pi/4; 
Phi_eta=-pi/2; 
% % Phi_eta=-pi/2;     % rad.����Phi_eta��[-pi,pi];����Pitching phase offset�����ȽϺ�
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % eta_m=0;                       % rad.��Pitching amplitude��Min
% % eta_m=pi/2;  
% % eta_m=pi;                     % rad.��Pitching amplitude��Max
% % %%%%%%%%%%%%%%
% % eta_0=eta_m-pi;            % rad.��Pitching offset��Min
% % eta_0=pi-eta_m;            % rad.��Pitching offset��Max
%%%%%%%%%%%%%%%%%%%%%%%%%
% K=0.0001;
% K=0.95; 
% K=1;       % K��(0,1];����K�ӽ�0ʱphiΪ��������; K=1ʱphiΪ���Ƿ���; %K�ɱ������ǳ������淴�䷽��Ķ��������淴�ٶȵĶ�����
%%%%%%%%%%%%%%%%%%%%%%%%%
% C_eta=inf; % C_eta��(0,inf];����C_eta�ӽ�0ʱetaΪ��������; C_eta=infʱeta����Ϊ��Ծ����; % C_eta��ֵ�������淴��ʱ���ʸ���أ�����
% C_eta=2.375; 
%%%%%%%%%%%%%%%%%%%%%%%%%
% C_eta=0.0001;
% C_eta=15; 
% C_eta=30; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% w =1185.6;     % ��Ƶ��   %  f=188.7; T=1/f;  %����Ƶ�� (Hz)������ 
% f=188.7; 
T=1/f;  %����Ƶ�� (Hz)������  % w =1185.6; 
% t_00=solve(dphi(phi_m,K,f,t)=0,t); 
t_00= 1/(4*f);    % t_00= -1/(4*f);
t=linspace(t_00,t_00+T,1000); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%7����������(f,phi_m,K,eta_m,C_eta,Phi_eta,eta_0)
K=[0.0001,0.95,1];                                              % Dickinson����ѡ����K_phi=0.01;
C_eta=[0.0001,5,10000];                                     % Dickinson����ѡ����C_eta=1.5;
n=length(t);
phi=zeros(n,3);
psi=zeros(n,4);
for i=1:length(K)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi(:,i)=phi_m.*asin(K(:,i).*sin(2.*pi.*f.*t))./asin(K(:,i));    % �Ĵ�ǡ����ٶȺͽǼ��ٶ�
% dphi=(2.*K.*pi.*f.*phi_m.*cos(2.*pi.*f.*t))./(asin(K).*(1 - K.^2.*sin(2.*f.*pi.*t).^2).^(1./2));
% ddphi=(4.*K.^3.*pi.^2.*f.^2.*phi_m.*cos(2.*pi.*f.*t).^2.*sin(2.*pi.*f.*t))./(asin(K).*(1 - K.^2.*sin(2.*f.*pi.*t).^2).^(3./2))...
%             - (4.*K.*pi.^2.*f.^2.*phi_m.*sin(2.*pi.*f.*t))./(asin(K).*(1 - K.^2.*sin(2.*f.*pi.*t).^2).^(1./2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psi(:,i)=-(eta_m.*tanh(C_eta(:,i).*sin(2.*pi.*f.*t+Phi_eta))./tanh(C_eta(:,i))+eta_0); % Ťת�ǡ����ٶȺͽǼ��ٶȡ���ԭʼ��ʽû�и���;
% dpsi=(2.*C_eta.*pi.*eta_m.*f.*cos(Phi_eta + 2.*pi.*f.*t).*(tanh(C_eta.*sin(Phi_eta + 2.*pi.*f.*t)).^2 - 1))./tanh(C_eta);% Ťת���ٶ�
% ddpsi=- (4.*C_eta.*pi.^2.*eta_m.*f.^2.*sin(Phi_eta + 2.*pi.*f.*t).*(tanh(C_eta.*sin(Phi_eta + 2.*pi.*f.*t)).^2 - 1))./tanh(C_eta)...
%       - (8.*C_eta.^2.*pi.^2.*eta_m.*f.^2.*cos(Phi_eta + 2.*pi.*f.*t).^2.*tanh(C_eta.*sin(Phi_eta + 2.*pi.*f.*t)).*(tanh(C_eta.*sin(Phi_eta + 2.*pi.*f.*t)).^2 - 1))./tanh(C_eta);% Ťת�Ǽ��ٶ�
end
C_eta=2.5;
Phi_eta=0; 
psi(:,4)=-(eta_m.*tanh(C_eta.*sin(2.*pi.*f.*t+Phi_eta))./tanh(C_eta)+eta_0);
% size(phi)
% size(psi)
%%%%%%%%%%%%%%%
[phimax,index1]=max(phi);
phi_max =phimax*180/pi;
% t_01=t(index1)
%%%%%%%%%%%%%%%
[psimax,index2]=max(psi);
psi_max =psimax*180/pi;
% t_02=t(index2)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)       % ͼ1�����Ĵ�Ǻ�Ťת��
hold on
subplot(211)
hold on
h1=plot(t/T,phi(:,1)*180/pi/phi_max(1,2),'r:',t/T,phi(:,2)*180/pi/phi_max(1,2),'r-.',t/T,phi(:,3)*180/pi/phi_max(1,2),'r-','LineWidth',2); %ת��Ϊms �� ����degree   *10^3   *180/pi
% xlabel('Normalized time')
ylabel('Normalized flapping angle','FontSize',14,'FontName','Times','FontWeight','Bold')
% legend('\it\phi(t)')
legend('\itK_\phi\rm= 0.0001','\itK_\phi\rm= 0.95','\itK_\phi\rm= 1')
% title('�Ĵ�Ǻ�Ťת����ʱ��ı仯����')   % �Ĵ�Ǻ�Ťת����ʱ��ı仯����
% grid on
axis([min(t)/T,max(t)/T,-1.2,1.2]) 
box on
set(gca,'LineStyle','-','LineWidth',1.5,'FontSize',12,'FontName','Times','FontWeight','Bold') 
% set(gca,'Color','r','LineStyle','-','LineWidth',3)

subplot(212)
hold on
h2=plot(t/T,psi(:,1)*180/pi/psi_max(1,2),'b:',t/T,psi(:,2)*180/pi/psi_max(1,2),'b-.',t/T,psi(:,3)*180/pi/psi_max(1,2),'b-','LineWidth',2); %ת��Ϊms �� ����degree   *10^3   *180/pi
xlabel('Normalized time','FontSize',18,'FontName','Times','FontWeight','Bold')
ylabel('Normalized pitch angle','FontSize',14,'FontName','Times','FontWeight','Bold')
% legend('\it\psi(t)')
hold on
plot(t/T,psi(:,4)*180/pi/psi_max(1,2),'g-','LineWidth',2)
legend('\itC_\psi\rm= 0.0001, \it\zeta\rm= -\pi/2','\itC_\psi\rm=5, \it\zeta\rm= -\pi/2','\itC_\psi\rm=10000, \it\zeta\rm= -\pi/2','\itC_\psi\rm=2.5, \it\zeta\rm= 0')
% title('�Ĵ�Ǻ�Ťת����ʱ��ı仯����')   % �Ĵ�Ǻ�Ťת����ʱ��ı仯����
% grid on
axis([min(t)/T,max(t)/T,-1.2,1.2]) 
box on
set(gca,'LineStyle','-','LineWidth',1.5,'FontSize',12,'FontName','Times','FontWeight','Bold') 
% set(gca,'Color','r','LineStyle','-','LineWidth',3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
