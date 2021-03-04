%读入数据
f1 = fopen('nT_4_nR_4_MMSE_FER.txt');
nT_4_nR_4_MMSE = textscan(f1,'%f %f %f %f %f %f');
f1 = fopen('nT_4_nR_4_V-BLAST_layers.txt');
nT_4_nR_4_V_BLAST_layers = textscan(f1,'%f %f %f %f %f %f %f %f %f');
f1 = fopen('nT_4_nR_6_ZF_MMSE.txt');
nT_4_nR_6_ZF_MMSE = textscan(f1,'%f %f %f %f %f %f %f %f %f %f');
%创建分方法测试横坐标
dB20 = 10.*log10(nT_4_nR_6_ZF_MMSE{1});
%创建V-BLAST分层测试横坐标
dB30 = 10.*log10(nT_4_nR_4_MMSE{1});
%作图
figure(1);
semilogy(dB30,nT_4_nR_4_MMSE{2},'--','LineWidth',2);grid on;hold on;
semilogy(dB30,nT_4_nR_4_MMSE{3},'o--','LineWidth',2);grid on;hold on;
semilogy(dB30,nT_4_nR_4_MMSE{4},'x-','LineWidth',2);grid on;hold on;
semilogy(dB30,nT_4_nR_4_MMSE{5},'x-.','LineWidth',2);grid on;hold on;
semilogy(dB30,nT_4_nR_4_MMSE{6},'x:','LineWidth',2);grid on;hold on;
title('n_T = 4 且 n_R = 4 的MMSE误帧率模拟实验，每次传输8Bits');
xlabel('$\frac{E_b}{N_0}$ in dB','Interpreter','latex');ylabel('FER');set(gca,'FontSize',22);
legend('MMSE-Linear(Pseudo inverse)','MMSE-V-BLAST','MMSE-Unsorted QRD','MMSE-SQRD','MMSE-SQRD-PSA');
figure(2);
semilogy(dB20,nT_4_nR_4_V_BLAST_layers{2},'s--','LineWidth',2);grid on;hold on;
semilogy(dB20,nT_4_nR_4_V_BLAST_layers{3},'o--','LineWidth',2);grid on;hold on;
semilogy(dB20,nT_4_nR_4_V_BLAST_layers{4},'x--','LineWidth',2);grid on;hold on;
semilogy(dB20,nT_4_nR_4_V_BLAST_layers{5},'d--','LineWidth',2);grid on;hold on;
semilogy(dB20,nT_4_nR_4_V_BLAST_layers{6},'s-','LineWidth',2);grid on;hold on;
semilogy(dB20,nT_4_nR_4_V_BLAST_layers{7},'o-','LineWidth',2);grid on;hold on;
semilogy(dB20,nT_4_nR_4_V_BLAST_layers{8},'x-','LineWidth',2);grid on;hold on;
semilogy(dB20,nT_4_nR_4_V_BLAST_layers{9},'d-','LineWidth',2);grid on;hold on;
title('n_T = 4 且 n_R = 4 的V-BLAST和MMSE-V-BLAST分层模拟实验，每次传输8Bits');
xlabel('$\frac{E_b}{N_0}$ in dB','Interpreter','latex');ylabel('BER');set(gca,'FontSize',22);
legend('ZF-Layer 1','ZF-Layer 2','ZF-Layer 3','ZF-Layer 4',...
    'MMSE-Layer 1','MMSE-Layer 2','MMSE-Layer 3','MMSE-Layer 4','Location','SW');
figure(3);
semilogy(dB20,nT_4_nR_6_ZF_MMSE{2},'s--','LineWidth',2);grid on;hold on;
semilogy(dB20,nT_4_nR_6_ZF_MMSE{3},'o--','LineWidth',2);grid on;hold on;
semilogy(dB20,nT_4_nR_6_ZF_MMSE{4},'x--','LineWidth',2);grid on;hold on;
semilogy(dB20,nT_4_nR_6_ZF_MMSE{5},'d--','LineWidth',2);grid on;hold on;
semilogy(dB20,nT_4_nR_6_ZF_MMSE{6},'s-','LineWidth',2);grid on;hold on;
semilogy(dB20,nT_4_nR_6_ZF_MMSE{7},'o-','LineWidth',2);grid on;hold on;
semilogy(dB20,nT_4_nR_6_ZF_MMSE{8},'x-','LineWidth',2);grid on;hold on;
semilogy(dB20,nT_4_nR_6_ZF_MMSE{9},'d-','LineWidth',2);grid on;hold on;
semilogy(dB20,nT_4_nR_6_ZF_MMSE{10},'*-','LineWidth',2);grid on;hold on;
title('n_T = 4 且 n_R = 6 的ZF和MMSE误码率模拟实验，每次传输8Bits');
xlabel('$\frac{E_b}{N_0}$ in dB','Interpreter','latex');ylabel('BER');set(gca,'FontSize',22);
legend('ZF-Pseudo inverse','ZF-V-BLAST','ZF-Unsorted QRD','ZF-SQRD',...
    'MMSE-Linear(Pseudo inverse)','MMSE-V-BLAST','MMSE-Unsorted QRD','MMSE-SQRD','MMSE-SQRD-PSA',...
    'Location','SW');
%关闭文件
fclose('all');