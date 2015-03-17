 H = 10;
 W = 10;
l1 = 5; 
l2 = 2;
% 
I1 = zeros(H,W);
I2 = zeros(H,W);
I1(2,3) = 200;
I1(2,2) = 200;
I1(3,3) = 200;
I1(3,2) = 200;

I2(4,3) = 200;
I2(4,2) = 200;
I2(3,3) = 200;
I2(3,2) = 200;


ix = conv2(I1,0.25* [-1 1; -1 1],'same') + conv2(I2, 0.25*[-1 1; -1 1],'same');
iy = conv2(I1, 0.25*[-1 -1; 1 1], 'same') + conv2(I2, 0.25*[-1 -1; 1 1], 'same');
ft = conv2(I1, 0.25*ones(2),'same') + conv2(I2, -0.25*ones(2),'same');
    



E0 = sparse(1:W,1:W,horzcat(-1/2,zeros(1,W-2),1/2),W,W);
E1 = sparse(2:W,1:W-1,horzcat([1/2,1/2],1/8*ones(1,W-4),1/2),W,W);
E2 = sparse(2:W,1:W-1,horzcat([1/2],1/8*ones(1,W-4),[1/2,1/2]),W,W);
E3 = sparse(3:W,1:W-2,horzcat([0,0],1/12*ones(1,W-4)),W,W);
E4 = sparse(3:W,1:W-2,horzcat(1/12*ones(1,W-4),[0,0]),W,W);

DW = E0+E1-E2'-E3+E4';


E0 = sparse(1:H,1:H,horzcat(1/2,zeros(1,H-2),-1/2),H,H);
E1 = sparse(2:H,1:H-1,horzcat(1/2,1/8*ones(1,H-4),1/2,1/2),H,H);
E2 = sparse(2:H,1:H-1,horzcat([1/2,1/2],1/8*ones(1,H-4),1/2),H,H);
E3 = sparse(3:H,1:H-2,horzcat([0,0],1/12*ones(1,H-4)),H,H);
E4 = sparse(3:H,1:H-2,horzcat(1/12*ones(1,H-4),[0,0]),H,H);

DH= E0+E1-E2'+E3'-E4;


A6D = zeros(H,W,l1,H,W,l2);

for i = 1:H
    for j = 1:W
        A6D(i,j,1,i,j,1) = ix(i,j);
        A6D(i,j,1,i,j,2) = iy(i,j);
        for k = 1:W
             A6D(i,j,2,i,k,1) = DW(k,j);
            A6D(i,j,3,i,k,2) = DW(k,j);
        end
        for k = 1:H
             A6D(i,j,4,k,j,1) = DH(k,i);
             A6D(i,j,5,k,j,2) = DH(k,i);
        end
    end
end

B = zeros(W,H,5);
B(:,:,1) = ft';
B0 = sparse(reshape(B,W*H*5,1));

A = sparse(reshape(A6D,[l1*H*W,l2*H*W]));

spy(A);
AOrt  = null(full(A)');
 
B = ft';
A1 = horzcat(AOrt',-1*AOrt');
bw = AOrt' * B0;
Asp = sparse(A1);
w0 = linprog([],[],[],Asp,bw,zeros(H*W*5),[]);

Im = eye(H*W*5);
bw0 = horzcat(Im,(-1)*Im) * w0 - B0;
x0 = mldivide(A, bw0);
x =reshape(x0, W, H, 2);
imwrite(uint8(flowToColor(x)), 'flow.bmp'); 

