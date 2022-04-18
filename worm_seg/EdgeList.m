function edgeList = EdgeList(data)
% ������������ȡ��ֵͼ�ı�Ե����
% dataΪ��ֵͼ��1ΪĿ�꣬0Ϊ����
% edgeList Ϊdata�ı�Ե����2*N�����������ĩ����Ϊͬһ��Ե��

direction = [0 1; -1 1; -1 0; -1 -1; 0 -1; 1 -1; 1 0; 1 1];
direction = [direction; direction];

% ��չ���ݣ�ʹ���ܱ�����Ϊ0
im = zeros(size(data,1) + 2, size(data,2) + 2);
im(2:(end-1), 2:(end-1)) = data;

%������ʼ��
[x, y] = find(im == 1);
xori = x(1); yori = y(1);

edgex = xori;
edgey = yori;

% ������ʼ���ҷ���
i = 6;
while 1
    if (1 == im(edgex(end) + direction(i, 1), edgey(end) + direction(i, 2)))
        edgex = [edgex, edgex(end) + direction(i, 1)];
        edgey = [edgey, edgey(end) + direction(i, 2)];

        offset = 2 - mod(i, 2);
        i = mod( (i - offset) - 1, 8) + 1 ;     
        if(edgex(1) == edgex(end) && edgey(1) == edgey(end))
            break;
        end
    else
        i = i+1;
    end    
end

edgeList = [edgex-1; edgey-1];
end
