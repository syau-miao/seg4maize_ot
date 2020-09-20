%本代码是论文《基于骨架和最优传输距离的玉米点云茎叶分割和表型提取》的测试代码
%我们只对论文中使用的30个样本进行了测试，如果您在运行时有其他的问题，可以联系
%代码的作者苗腾，邮箱为miaoteng@syau.edu.cn
%我们的代码只是一种技术思路的简单实现，没有做过优化和严格的测试，因此不保证代码
%在处理其他数据时会没有bug。
%开展农业工程研究的各位同仁可以随意修改和使用该代码，本人使用Matlab时间较短，水平较差，
%代码有不堪入目的地方请见谅。
%特别感谢大连理工大学的曹俊杰老师，在github上分享了骨架提取的代码，我们的骨架
%提取采用了曹老师的这部分工作。

clear;clc;close all;
path('toolbox',path);
path('skeleton',path);
path('PhenotypicTrait',path);
path('segment',path);
path('Evaluate',path);
warning('off');
tic;
IdName='30';%我们共有30个数据，文件名分别为1到30
            %我们的数据都是已经分割好的文件，每个
            %器官一个Ply,一个植株在一个文件夹下。
            %%在进行分割时会合到一起。
filename=['testdata\' IdName '\'];
DebugShow=true;
%************ 用户参数初始化***********************************
Parameters.KnnNum=16; %%%% 论文中用来控制骨架局部细节的K值 ( paper 2.3.1 )
Parameters.t2=0.2;
Parameters.alpha=1.0;%%%论文中的alpha 用于控制茎点云的粗细和长度
isWrite=0;%%%%% 这个参数设为1 时，将根据用户在上面的设定执行算法
          %%%%% 这个参数设为0时，K参数和alpha将调用我们对每个点云文件设定的比较好的参数
          %%%%%% 但e参数仍由用户指定，便于我们更加直观地理解论文中e参数的作用
isWritePhe=0;
Parameters=SerializeParameter(Parameters,IdName,isWrite);
Parameters.e=5; %%%%%%%%论文中公式（2）的e，此处修改参数e
%************1 骨架提取*****************************************
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
%************* 2 基于骨架的粗分割 *********************
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
%********************3 精细分割 ***********************************
EMD=computeEMD(PA_Pts,5);
AutoRegions=classifyLeafPoints_EMD(PA_Pts,EMD,Phi_O,Phi_U,1,1);
%**************************************************************************
%%%%%%%%%%%%%%%%%%%表型提取%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
traitfile=[IdName '.txt'];%%%%表型提取的结果将存在一个txt文件中
PhenotypicTrait(PA_Pts,AutoRegions,traitfile,DebugShow);

