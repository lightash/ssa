function [sig] = testSignal(N, mAmp, mPh)
        
    if nargin < 2   % by default or without 2nd arg
        
        mas1 = 1;       % ��������� ��������
        mas2 = 0;%.5;
        mas3 = 0;%.25;
        mapi = 0;%.75;    % ��������� ����
        mame = 0;%.35;    %           �������
        a31 = 1;%.1;      % ������������� ���� �����
        a32 = 0;%.3;
        apost = 0;      % ���������� ������������
        tet1 = 0;       % ��������� ���� ��������
        tet2 = 0;
        tet3 = 0;
        
        mAmp = [mas1,mas2,mas3,mapi,mame,a31,a32,apost];
        mPh = [tet1,tet2,tet3];
        
    end
    
    t = 0:1/N:1-1/N;
    
    % ���������
    as1 = mAmp(1)*sin(2*pi*t+mPh(1));
    as2 = mAmp(2)*sin(2*pi*2*t+mPh(2));
    as3 = mAmp(3)*sin(2*pi*4*t+mPh(3));
    % ����
    api = mAmp(4)*( 2*pulstran(t,0:1/16:1,'tripuls',1/16,1)-1 );
    % ������
    ame(1:N) = mAmp(5);
    for i = 3:4:N
       ame(i) = -mAmp(5);
       ame(i+1) = -mAmp(5);
    end
    % ������������� �����
    a3(1:2:N) = mAmp(6);
    a3(2:2:N) = mAmp(7);

    sig = as1 + as2 + as3 + api + ame + a3 + mAmp(8);
    
end