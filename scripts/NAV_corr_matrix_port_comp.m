clc
% close all
clear all

load('D:\Dropbox\Signals\incartdb\I20\I20proc.mat')
in = val(1,:);
annot(2461) = 'V';   % Fusion of ventricular and normal beat

all_beats = 1:length(annot);
% Normal beat, Atrial premature beat, Premature ventricular contraction
bmark = 'NAV';
btypeN = 3;       % Number of beat types to examine
for btype = 1:btypeN
   Bnum{btype} = all_beats(annot == bmark(btype));                               % anniNAV
   Blen(btype) = length(Bnum{btype});
   Bpos{btype} = mark(Bnum{btype});                                              % ann
   Bord(all_beats(annot ==  bmark(btype))) = btype*ones(1,length(Bpos{btype}));  % annNAV
end
Bwin = [-47 80];   % Borders of PQRST period
winL = Bwin(2)-Bwin(1)+1;

% Generating portraits
perN = all_beats(end);  % Number of periods to use
f = zeros(perN, winL);
for per = 1:perN

   period = mark(per);
   window = period+Bwin(1): period+Bwin(2);
   f(per,:) = in(window);
   fm(per) = mean(f(per,:));
   f(per,:) = f(per,:) - fm(per);

   [f(per,:), fs(per)] = nrm(f(per,:));
end
port = cell(1,btypeN);
for btype = 1:btypeN
   port{btype} = nrm(mean(f(Bnum{btype},:),1));
end

btypes = {[1 2],[1 3],[2 3]};

%% Additional portraits
% % Добавление к портретам отсчётов среднего и масштаба для трёх пар
% f_aug = [f fm' fs'];
% for btype = 1:btypeN
%    port_aug{btype} = [port{btype} mean(fm) mean(fs)];
% end
% winL = 130;
% f = f_aug;
% port = port_aug;

% % Корреляционная матрица как портрет
% f_cormat = zeros(perN, winL^2/16);
% for per = 1:perN
%    f_cormat(per,:) = nrm(transform( f(per,1:4:end)' * f(per,1:4:end) ),1);
% end
% for btype = 1:btypeN
%    port_cormat{btype} = nrm(transform( port{btype}(1:4:end)' * port{btype}(1:4:end) ),1);
% end
% winL = winL^2/16;
% f = f_cormat;
% port = port_cormat;

% % Корреляционная функция как портрет
% f_xcorr = zeros(perN, 2*winL-1);
% for per = 1:perN
%    f_xcorr(per,:) = nrm(xcorr(f(per,:)),1);
% end
% for btype = 1:btypeN
%    port_xcorr{btype} = nrm(xcorr(port{btype}),1);
% end
% winL = 2*winL-1;
% f = f_xcorr;
% port = port_xcorr;

% % Суммарная корреляционная функция как портрет
% for per = 1:perN
%    cor = xcorr(f(per,:));
%    cor = cor(ceil(end/2):end);
%    cor(2:end) = 2*cor(2:end);
%    f_sumcor(per,:) = nrm(cor,1);
% end
% for btype = 1:btypeN
%    cor = xcorr(port{btype});
%    cor = cor(ceil(end/2):end);
%    cor(2:end) = 2*cor(2:end);
%    port_sumcor{btype} = nrm(cor,1);
% end
% f = f_sumcor;
% port = port_sumcor;

% % Комплексирование: портрет + портрет
% f_ff = zeros(perN, 2*winL);
% for per = 1:perN
%    f_ff(per,:) = nrm([ f(per,:) f(per,:) ]);
% end
% for btype = 1:btypeN
%    port_ff{btype} = nrm([ port{btype} port{btype} ]);
% end
% winL = 2*winL;
% f = f_ff;
% port = port_ff;

% % Комплексирование: портрет + корреляционная матрица
% f_fcormat = zeros(perN, winL + winL^2/16);
% for per = 1:perN
%    f_fcormat(per,:) = nrm([ nrm(transform( f(per,1:4:end)' * f(per,1:4:end) ),1) f(per,:) ],1);
% end
% for btype = 1:btypeN
%    port_fcormat{btype} = nrm([ nrm(transform( port{btype}(1:4:end)' * port{btype}(1:4:end) ),1) port{btype} ],1);
% end
% winL = winL + winL^2/16;
% f = f_fcormat;
% port = port_fcormat;

% Комплексирование: корреляционная функция + корреляционная матрица
coef = [sqrt(2*winL-1) sqrt(winL^2/16)]/(2*winL-1 + winL^2/16);
f_xcorrcormat = zeros(perN, 2*winL-1 + winL^2/16);
for per = 1:perN
   f_xcorrcormat(per,:) = ([ nrm(xcorr(f(per,:)),1)/coef(1) nrm(transform( f(per,1:4:end)' * f(per,1:4:end) ),1)/coef(2) ]);
end
for btype = 1:btypeN
   port_xcorrcormat{btype} = ([ nrm(xcorr(port{btype}),1)/coef(1) nrm(transform( port{btype}(1:4:end)' * port{btype}(1:4:end) ),1)/coef(2) ]);
end
winL = 2*winL-1 + winL^2/16;
f = f_xcorrcormat;
port = port_xcorrcormat;

% % Комплексирование: портрет + корреляционная функция
% coef = [sqrt(winL) sqrt(2*winL-1)]/(winL + 2*winL-1);
% f_xcorr = zeros(perN, winL + 2*winL-1);
% for per = 1:perN
%    f_xcorr(per,:) = ([f(per,:)/coef(1) nrm(xcorr(f(per,:)),1)/coef(2)]);
% end
% for btype = 1:btypeN
%    port_xcorr{btype} = ([port{btype}/coef(1) nrm(xcorr(port{btype}),1)/coef(2)]);
% end
% winL = winL + 2*winL-1;
% f = f_xcorr;
% port = port_xcorr;

% % Комплексирование: суммарная корреляционная функция как портрет
% for per = 1:perN
%    cor = xcorr(f(per,:));
%    cor = cor(ceil(end/2):end);
%    cor(2:end) = 2*cor(2:end);
%    cor = nrm(cor,1);
%    f_sumcor(per,:) = nrm([f(per,:) cor], 1);
% end
% for btype = 1:btypeN
%    cor = xcorr(port{btype});
%    cor = cor(ceil(end/2):end);
%    cor(2:end) = 2*cor(2:end);
%    cor = nrm(cor,1);
%    port_sumcor{btype} = nrm([port{btype} cor], 1);
% end
% winL = 2*winL;
% f = f_sumcor;
% port = port_sumcor;


%% Итеративный корреляционный отбор в трёх парах N, A, V с вероятностями и локализациями
for bpair = 1:3
   beats = sort([Bnum{btypes{bpair}(1)} Bnum{btypes{bpair}(2)}]);
   beatsN = Blen(btypes{bpair}(1))+Blen(btypes{bpair}(2));

   % Full probability
   des = zeros(btypeN);
   for per = 1:beatsN
      cor = -ones(1,3);
      cor(btypes{bpair}(1)) = f(beats(per),:)*port{btypes{bpair}(1)}';
      cor(btypes{bpair}(2)) = f(beats(per),:)*port{btypes{bpair}(2)}';
      [~,ind] = max(cor);
      des(Bord(beats(per)),ind) = des(Bord(beats(per)),ind) + 1/Blen(Bord(beats(per)));
   end
   prob_all{bpair}(:,1) = [des(btypes{bpair}(1),btypes{bpair}(1));des(btypes{bpair}(2),btypes{bpair}(2))];
   
   % Formal criteriums
   H = CritHist({ f(Bnum{btypes{bpair}(1)},:), f(Bnum{btypes{bpair}(2)},:) });
   I = DK( H{1}, H{2}, Blen(btypes{bpair}(1)), Blen(btypes{bpair}(2)) );
%    I = AlphaZ( H{1}, H{2} );
   [~,order] = sort(I,'descend');  % keep first informative ones

   % Excuding serially
   win_rem{bpair,1} = 1:winL;
   for throw = 1:winL-1
      disp([bpair throw])
      
%       % Energy criterium
%       prob_one = zeros(2,winL-throw+1);
%       for exclude = 1:winL-throw+1
% 
%          win = [win_rem{bpair,throw}(1:exclude-1) win_rem{bpair,throw}(exclude+1:end)];
% 
%          des = zeros(btypeN);
%          for per = 1:beatsN
%             cor = -ones(1,3);
%             cor(btypes{bpair}(1)) = nrm(f(beats(per),win),1) * nrm(port{btypes{bpair}(1)}(win),1)';
%             cor(btypes{bpair}(2)) = nrm(f(beats(per),win),1) * nrm(port{btypes{bpair}(2)}(win),1)';
%             [~,ind] = max(cor);
%             des(Bord(beats(per)),ind) = des(Bord(beats(per)),ind) + 1/Blen(Bord(beats(per)));
%          end
%          prob_one(:,exclude) = [des(btypes{bpair}(1),btypes{bpair}(1));des(btypes{bpair}(2),btypes{bpair}(2))];
%       end
%       [~,ind] = max(mean(prob_one,1));
%       prob_all{bpair}(:,throw+1) = prob_one(:,ind);
%       order{bpair}(throw) = find( port{btypes{bpair}(1)} == port{btypes{bpair}(1)}(win_rem{bpair,throw}(ind)), 1, 'first');
%       win_rem{bpair,throw+1} = [win_rem{bpair,throw}(1:ind-1) win_rem{bpair,throw}(ind+1:end)];
      
      % Formal criteriums
      win = sort(order(1:winL-throw));
      des = zeros(btypeN);
      for per = 1:beatsN
         cor = -ones(1,3);
         cor(btypes{bpair}(1)) = nrm(f(beats(per),win),1) * nrm(port{btypes{bpair}(1)}(win),1)';
         cor(btypes{bpair}(2)) = nrm(f(beats(per),win),1) * nrm(port{btypes{bpair}(2)}(win),1)';
         [~,ind] = max(cor);
         des(Bord(beats(per)),ind) = des(Bord(beats(per)),ind) + 1/Blen(Bord(beats(per)));
      end
      win_rem{bpair,throw+1} = win;
      prob_all{bpair}(:,throw+1) = [ des(btypes{bpair}(1),btypes{bpair}(1)); des(btypes{bpair}(2),btypes{bpair}(2)) ];
      
   end
   
%    % Energy criterium
%    order{bpair}(winL) = win_rem{bpair,end};
end

save('NAVcmpc_ku_xcorr+cormat_coef.mat','prob_all','win_rem','order')
%% Graphs
figure('color','white')
for bpair = 1:3
   [m,i] = max(mean(prob_all{bpair}(:,end:-1:1),1));
   i = winL - i + 1;
   
   subplot(1,3,bpair)
   plot(prob_all{bpair}(1,:),'k','LineWidth',1.5),hold on
   plot(prob_all{bpair}(2,:),'--k','LineWidth',1.5)
   plot(mean(prob_all{bpair},1),'-.k','LineWidth',1.5),grid,axis tight,ylim([.8 1])
   plot([i i],[.5 1],'--k','LineWidth',3)
   
   legend(bmark(btypes{bpair}(1)),bmark(btypes{bpair}(2)),'Среднее')
   title([m i])
end
subplot(131)
ylim([.7 .9+1e-9])
ylabel('Вероятности распознавания','FontName','Times New Roman','FontSize',14)
subplot(132)
xlabel('Количество исключённых малоинформативных отсчётов','FontName','Times New Roman','FontSize',14)

%% Localizations graphs
% x = (Bwin(1):Bwin(2))/Fd*1e3;
figure('color','white')
for bpair = 1:3
   [m,i] = max(mean(prob_all{bpair}(:,end:-1:1),1));
   i = winL - i + 1;
   
   subplot(1,3,bpair)
   plot(port{btypes{bpair}(1)},'-k'),hold on%x,
   plot(port{btypes{bpair}(2)},'--k')%x,
   plot(win_rem{bpair,i},port{btypes{bpair}(1)}(win_rem{bpair,i}),'.k','LineWidth',1.5)%/Fd*1e3 + x(1)-1e3/Fd
   plot(win_rem{bpair,i},port{btypes{bpair}(2)}(win_rem{bpair,i}),'.k','LineWidth',1.5),grid,axis tight%/Fd*1e3 + x(1)-1e3/Fd
   
   legend(bmark(btypes{bpair}(1)),bmark(btypes{bpair}(2)))
   title([m i])
end
subplot(131)
ylabel({'Характеристики формы' 'QRS-комплексов по типам'},'FontName','Times New Roman','FontSize',14)
subplot(132)
% xlabel('Время от R-зубца, мс','FontName','Times New Roman','FontSize',14)
xlabel('Номер отсчёта','FontName','Times New Roman','FontSize',14)

% %% Гистограммы корр. интеграла функций с форм. критериями
% for bpair = 1:3
%    beats = sort([Bnum{btypes{bpair}(1)} Bnum{btypes{bpair}(2)}]);
%    beatsN = Blen(btypes{bpair}(1))+Blen(btypes{bpair}(2));
% 
%    % Full probability
%    des = zeros(btypeN);
%    corrs{bpair} = -ones(beatsN,3);
%    for per = 1:beatsN
%       corrs{bpair}(per,btypes{bpair}(1)) = f(beats(per),:)*port{btypes{bpair}(1)}';
%       corrs{bpair}(per,btypes{bpair}(2)) = f(beats(per),:)*port{btypes{bpair}(2)}';
%       [~,ind] = max(corrs{bpair}(per,:));
%       des(Bord(beats(per)),ind) = des(Bord(beats(per)),ind) + 1/Blen(Bord(beats(per)));
%    end
%    
%    H = CritHist({ corrs{bpair}(:,btypes{bpair}(1)), corrs{bpair}(:,btypes{bpair}(2)) });
%    Ik(bpair) = DK( H{1}, H{2}, Blen(btypes{bpair}(1)), Blen(btypes{bpair}(2)) );
%    Ia(bpair) = AlphaZ( H{1}, H{2} );
% end
% 
% figure('color','white')
% subplot(211),stem(Ik,'k','LineWidth',1.5),title('Kullback'),grid,xlim([0 4])
% subplot(212),stem(Ia,'k','LineWidth',1.5),title('\alpha_z'),grid,xlim([0 4])

% %% Гистограммы пар "корр. матр, корр. функции (ро) и сумм. корр. ф. (эр)" с форм. критериями
% % Корреляционная матрица как портрет
% fcor{1} = zeros(perN, winL^2);
% for per = 1:perN
%    fcor{1}(per,:) = nrm(transform( f(per,:)' * f(per,:) ),1);
% end
% 
% % Корреляционная функция как портрет
% fcor{2} = zeros(perN, 2*winL-1);
% for per = 1:perN
%    fcor{2}(per,:) = nrm(xcorr(f(per,:)),1);
% end
% 
% % Суммарная корреляционная функция как портрет
% fcor{3} = zeros(perN, winL);
% for per = 1:perN
%    cor = xcorr(f(per,:));
%    cor = cor(ceil(end/2):end);
%    cor(2:end) = 2*cor(2:end);
%    fcor{3}(per,:) = nrm(cor,1);
% end
% 
% for pair = 1:3
%    H = CritHist({ fcor{pair}(:,btypes{pair}(1)), fcor{pair}(:,btypes{pair}(2)) });
%    Ik(pair) = DK( H{1}, H{2}, Blen(btypes{pair}(1)), Blen(btypes{pair}(2)) );
%    Ia(pair) = AlphaZ( H{1}, H{2} );
% end
% 
% figure('color','white')
% subplot(211),stem(Ik,'k','LineWidth',1.5),title('Kullback'),grid,xlim([0 4])
% subplot(212),stem(Ia,'k','LineWidth',1.5),title('\alpha_z'),grid,xlim([0 4])
% 

















