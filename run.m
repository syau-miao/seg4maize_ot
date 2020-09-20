%�����������ġ����ڹǼܺ����Ŵ����������׵��ƾ�Ҷ�ָ�ͱ�����ȡ���Ĳ��Դ���
%����ֻ��������ʹ�õ�30�����������˲��ԣ������������ʱ�����������⣬������ϵ
%������������ڣ�����Ϊmiaoteng@syau.edu.cn
%���ǵĴ���ֻ��һ�ּ���˼·�ļ�ʵ�֣�û�������Ż����ϸ�Ĳ��ԣ���˲���֤����
%�ڴ�����������ʱ��û��bug��
%��չũҵ�����о��ĸ�λͬ�ʿ��������޸ĺ�ʹ�øô��룬����ʹ��Matlabʱ��϶̣�ˮƽ�ϲ
%�����в�����Ŀ�ĵط�����¡�
%�ر��л��������ѧ�Ĳܿ�����ʦ����github�Ϸ����˹Ǽ���ȡ�Ĵ��룬���ǵĹǼ�
%��ȡ�����˲���ʦ���ⲿ�ֹ�����

clear;clc;close all;
path('toolbox',path);
path('skeleton',path);
path('PhenotypicTrait',path);
path('segment',path);
path('Evaluate',path);
warning('off');
tic;
IdName='30';%���ǹ���30�����ݣ��ļ����ֱ�Ϊ1��30
            %���ǵ����ݶ����Ѿ��ָ�õ��ļ���ÿ��
            %����һ��Ply,һ��ֲ����һ���ļ����¡�
            %%�ڽ��зָ�ʱ��ϵ�һ��
filename=['testdata\' IdName '\'];
DebugShow=true;
%************ �û�������ʼ��***********************************
Parameters.KnnNum=16; %%%% �������������ƹǼֲܾ�ϸ�ڵ�Kֵ ( paper 2.3.1 )
Parameters.t2=0.2;
Parameters.alpha=1.0;%%%�����е�alpha ���ڿ��ƾ����ƵĴ�ϸ�ͳ���
isWrite=0;%%%%% ���������Ϊ1 ʱ���������û���������趨ִ���㷨
          %%%%% ���������Ϊ0ʱ��K������alpha���������Ƕ�ÿ�������ļ��趨�ıȽϺõĲ���
          %%%%%% ��e���������û�ָ�����������Ǹ���ֱ�۵����������e����������
isWritePhe=0;
Parameters=SerializeParameter(Parameters,IdName,isWrite);
Parameters.e=5; %%%%%%%%�����й�ʽ��2����e���˴��޸Ĳ���e
%************1 �Ǽ���ȡ*****************************************
%%%% default parameter from paper(Cao et al 2010) 
sampleScale=0.02;
t1 = 0.1; % for inner branch nodes
a1 = pi*5.0/7.0; % for inner branch nodes, 
t2 = Parameters.t2;    
t3 = 30; % for small cycles;
[P.pts,trueRegions]=loadSegmentFile2(filename); 
eg_skeleton_laplacian_rosa;
eg_refine_skeleton;
%**************************************************************************
%************* 2 ���ڹǼܵĴַָ� *********************
[joints ,roots, branches]=find_Joints_mt(P.spls, P.spls_adj,P.pts, false); 
[sub_skeletons,P.spls]=SkeletonDecomposition(P.spls, joints,roots, P.spls_adj,false);
[Phi_O,Phi_U]=CoaseSegBySkeleton(P.pts,P.corresp,sub_skeletons,joints,3,false);
%%%% transform global coordinate axis to plant coordinate axis
[PA_Pts,PA_Spls]= ConstructPlantAxis(P.pts,P.spls,Phi_U,sub_skeletons);
if(DebugShow)
    ShowSkeleton4Paper(PA_Spls,P.spls_adj,joints ,roots, branches,PA_Pts,Phi_O,Phi_U);
end
%%%%% stem constraint 
[PA_StopX,PA_StartX]= FindLeafStopZ(sub_skeletons,joints,PA_Spls);  
[Phi_S,Phi_U]=StemSegment(PA_Pts,Phi_U,PA_StartX,Parameters.alpha,DebugShow);
Phi_O=flip(Phi_O);
Phi_O{end+1}=Phi_S;
PA_StopX(end+1)=-inf;
PA_StopX=flip(PA_StopX);
Phi_O=flip(Phi_O);
%**************************************************************************
%********************3 ��ϸ�ָ� ***********************************
EMD=computeEMD(PA_Pts,5);
AutoRegions=classifyLeafPoints_EMD(PA_Pts,EMD,Phi_O,Phi_U,1,1);
%**************************************************************************
%%%%%%%%%%%%%%%%%%%������ȡ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
traitfile=[IdName '.txt'];%%%%������ȡ�Ľ��������һ��txt�ļ���
PhenotypicTrait(PA_Pts,AutoRegions,traitfile,DebugShow);

