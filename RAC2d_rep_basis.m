function [monomials NewR reps] = RAC2d_rep_basis(d)
% Constructs the irreducible decomposition
%
% Based on original code from Marc-Olivier Renou
%
% The representations appear in the same order as in Table IV of
% the companion paper. Refer to that table for the dimensions and
% multiplicities.
%
% The final multiplicities are:
%
% - 5 representations T
% - 3 representations S
% - 7 representations Phi
% - 4 representations pi+
% - 3 representations pi-
% - 1 representations Lambda
% - 1 representations Omega
% - 1 representations lambda
% - 1 representations omega

    assert(d > 2, 'Supports only d > 2');

    monomials = cell(1, 0);
    monomials{1} = [];

    count = 1;
    for x1 = 1:d
        for x2 = 1:d
            count = count + 1;
            monomials{count} = x1 + (x2 - 1)*d;
        end
    end

    for y = 1:2 % n=2
        for b = 1:d
            count = count + 1;
            monomials{count} = d^2 + b + (y - 1)*d;
        end
    end


    for x1 = 1:d
        for x2 = 1:d
            for y = 1:2 %n=2
                for b = 1:d
                    count = count + 1;
                    monomials{count} = [x1+(x2-1)*d d^2+b+(y-1)*d];
                end
            end
        end
    end

    reps = [1 1 2*(d-1) (d-1)^2 (d-1)^2 (d-1)*(d-2) d*(d-3) (d-1)^2*(d-2) d*(d-1)*(d-3)
            5 3 7       4       3       1           1       1             1];
    reps = reps(:, (reps(1,:)>0) & (reps(2,:)>0));
    n=2;

    G=[];

    MatrixDim=1+d^2+2*d+2*d^3;
    G=zeros(MatrixDim^2,2);

    id=eye(d);
    %SymOperators=SymmetryConstructor(n,d); % Ug (group elements)


    %% Block diagonal basis creation


    deltaY(:,1)=[1 1]'/sqrt(2);
    deltaY(:,2)=[1 -1]'/sqrt(2);
    id2=eye(2);

    b0=ones(1,d)'/sqrt(d);

    for k=1:d-1
        dummy=0;
        for ii=1:k
            dummy=dummy+id(:,ii);
        end
        dummy(k+1)=-k;
        xi(:,k)=dummy/norm(dummy);
    end

    xi=[b0 xi]; % Basis elements



    for i=1:d
        for j=1:d
            delta(:,i,j)=id(:,i)-id(:,j);
        end
    end



    %%%%%%%%%%%% Basis for W (antisymmetric square of standard representation of the symmetric group) %%%%%%%%%%%%%%%%%

    Wbeta=[];
    loc=0;
    for j=2:d
        for i=2:j-1
            loc=loc+1;
            Wbeta=[Wbeta (Tensor(xi(:,i),xi(:,j))-Tensor(xi(:,j),xi(:,i)))/norm(Tensor(xi(:,i),xi(:,j))-Tensor(xi(:,j),xi(:,i)))];
        end
    end
    Wbeta;


    %%%%%%%%%%%%%%%%%%% Basis for T (trivial part of the symmetric square of the standard rep of the sym group) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Tt=0;
    for j=1:d
        for i=1:j-1
            Tt=Tt+Tensor(delta(:,i,j),delta(:,i,j));
        end
    end
    Tt=Tt/norm(Tt);

    %%%%%%%%%%%%% Basis for Delta (standard rep part of the sym square of the standard rep of the sym group) %%%%%%%%%%%%%%%%%%%%%%

    for k=1:d
        dummy=0;
        for i=1:d
            if i==k
            else
                dummy=dummy+Tensor(delta(:,i,k),delta(:,i,k));
            end
        end
        
        for j=1:d
            for i=1:j-1
                if i==k || j==k
                    [i j k];
                else
                    dummy=dummy-Tensor(delta(:,i,j),delta(:,i,j));
                end
            end
        end
        
        dummy=dummy+(d-4)/d*Tt;
        
        standardDelta(:,k)=dummy;
        
    end



    for k=1:d
        canonicalbasis(:,k)=standardDelta(:,k)+sqrt(8)*Tt;
        canonicalbasis(:,k)=canonicalbasis(:,k)/norm(canonicalbasis(:,k));
    end
    for k=1:d-1
        dummy=0;
        for ii=1:k
            dummy=dummy+canonicalbasis(:,ii);
        end
        dummy=dummy-k*canonicalbasis(:,k+1);
        Dxi(:,k)=dummy/norm(dummy);
    end



    %%%%%%%%%%% Psi %%%%%%%%%%%%%%%%%%%

    DeltaDeltaid=0;
    for i=2:d
        for j=2:d
            DeltaDeltaid=DeltaDeltaid+Tensor(xi(:,i),xi(:,j))*Tensor(xi(:,i),xi(:,j))';
        end
    end


    PsiPsiid=DeltaDeltaid;
    for i=1:size(Wbeta,2)
        PsiPsiid=PsiPsiid-Wbeta(:,i)*Wbeta(:,i)';
    end
    PsiPsiid=PsiPsiid-Tt*Tt';
    for i=1:size(Dxi,2)
        PsiPsiid=PsiPsiid-Dxi(:,i)*Dxi(:,i)';
    end

    Psiu = [];
    if d > 3
        Psiu = orth(PsiPsiid);
    end


    O=[Wbeta Tt Dxi Psiu];



    %%%%%%%%%%%%%%%%%%%%%%% X part
    R1=[];
    XbasisT=Tensor(xi(:,1),xi(:,1)); % trivial representation
    R1(:,1)=XbasisT;

    loccount=1;
    for k=2:d
        loccount=loccount+1;
        dummy=Tensor(xi(:,1),xi(:,k)); % phi rep
        XbasisPhi(:,k-1,1)=dummy;
        R1(:,loccount)=dummy;
    end

    for k=2:d
        loccount=loccount+1;
        dummy=Tensor(xi(:,k),xi(:,1)); % phi rep
        XbasisPhi(:,k-1,2)=dummy;
        R1(:,loccount)=dummy;
    end

    for k=2:d
        for l=2:d
            XbasisS3(:,k-1,l-1)=Tensor(xi(:,k),xi(:,l)); % s3 rep
            loccount=loccount+1;
            R1(:,loccount)=XbasisS3(:,k-1,l-1);
        end
    end

    
    %%%%%%%%%%%%%%%%%% yb part
    R2=[];
    YBbasisT=Tensor(deltaY(:,1),xi(:,1)); % trivial rep
    YBbasisS2=Tensor(deltaY(:,2),xi(:,1)); % s2 rep
    R2(:,1)=YBbasisT;
    R2(:,2)=YBbasisS2;

    loccount=2;

    for k=2:d
        loccount=loccount+1;
        YBbasisPhi(:,k-1,2)=Tensor(id2(:,2),xi(:,k)); % phi rep
        R2(:,loccount)=YBbasisPhi(:,k-1,2);
    end

    for k=2:d
        loccount=loccount+1;
        YBbasisPhi(:,k-1,1)=Tensor(id2(:,1),xi(:,k)); % phi rep
        R2(:,loccount)=YBbasisPhi(:,k-1,1);
    end




    %%%%%%%%%%%%%%%%% AB part

    %% All exept S3Phi and PhiPhi

    R3=[];
    ABbasisTT=Tensor(XbasisT,YBbasisT);
    ABbasisTS2=Tensor(XbasisT,YBbasisS2);
    R3(:,1)=ABbasisTT;
    R3(:,2)=ABbasisTS2;

    loccount=2;
    for l=2:-1:1
        for k=1:d-1
            loccount=loccount+1;
            ABbasisTPhi(:,k,l)=Tensor(XbasisT,YBbasisPhi(:,k,l));
            R3(:,loccount)=ABbasisTPhi(:,k,l);
        end
    end

    for l=1:2
        for k=1:d-1
            loccount=loccount+1;
            ABbasisPhiT(:,k,l)=Tensor(XbasisPhi(:,k,l),YBbasisT);
            R3(:,loccount)=ABbasisPhiT(:,k,l);
        end
    end

    for k=1:d-1
        for l=1:d-1
            loccount=loccount+1;
            ABbasisS3T(:,k,l)=Tensor(XbasisS3(:,k,l),YBbasisT);
            R3(:,loccount)=ABbasisS3T(:,k,l);
        end
    end


    for k=1:d-1
        for l=1:d-1
            loccount=loccount+1;
            ABbasisS3S2(:,k,l)=Tensor(XbasisS3(:,k,l),YBbasisS2);
            R3(:,loccount)=ABbasisS3S2(:,k,l);
        end
    end

    for l=1:2
        for k=1:d-1
            loccount=loccount+1;
            ABbasisPhiS2(:,k,l)=(-1)^(l+1)*Tensor(XbasisPhi(:,k,l),YBbasisS2);%zeros(size(ABbasisS3S2(:,1,1)))%
            R3(:,loccount)=ABbasisPhiS2(:,k,l);
        end
    end




    %% S3Phi

    %S3Phi_W
    for k=1:size(Wbeta,2)
        for j=2:d
            r=1;
            dummy=zeros(2*d^3,1);
            for i=2:d
                for l=2:d
                    dummy=dummy+Tensor(xi(:,i),xi(:,l))'*Wbeta(:,k)*Tensor(xi(:,i),xi(:,j),id2(:,r),xi(:,l));
                end
            end
            loccount=loccount+1;
            ABbasisS3Phi_W(:,k,j,r)=dummy;
            R3(:,loccount)= ABbasisS3Phi_W(:,k,j,r);
        end
    end
    for k=1:size(Wbeta,2)
        for i=2:d
            r=2;
            dummy=zeros(2*d^3,1);
            for j=2:d
                for l=2:d
                    dummy=dummy+Tensor(xi(:,j),xi(:,l))'*Wbeta(:,k)*Tensor(xi(:,i),xi(:,j),id2(:,r),xi(:,l));
                end
            end
            loccount=loccount+1;
            ABbasisS3Phi_W(:,k,i,r)=dummy;
            R3(:,loccount)= ABbasisS3Phi_W(:,k,i,r);
        end
    end


    %S3Phi_t (equiv to Phi)
    for j=2:d
        r=1;
        dummy=zeros(2*d^3,1);
        for i=2:d
            for k=2:d
                dummy=dummy+Tensor(xi(:,i),xi(:,k))'*Tt*Tensor(xi(:,i),xi(:,j),id2(:,r),xi(:,k));
            end
        end
        loccount=loccount+1;
        ABbasisS3Phi_t(:,j,r)=dummy;
        R3(:,loccount)= ABbasisS3Phi_t(:,j,r);
        
    end
    for i=2:d
        r=2;
        dummy=zeros(2*d^3,1);
        for j=2:d
            for k=2:d
                dummy=dummy+Tensor(xi(:,j),xi(:,k))'*Tt*Tensor(xi(:,i),xi(:,j),id2(:,r),xi(:,k));
            end
        end
        loccount=loccount+1;
        ABbasisS3Phi_t(:,i,r)=dummy;
        R3(:,loccount)= ABbasisS3Phi_t(:,i,r);
    end



    %S3Phi_Delta


    r=1;
    for k=1:size(Dxi,2)
        for j=1:d-1
            dummy=zeros(2*d^3,1);
            for i=1:d-1
                for l=1:d-1
                    dummy=dummy+Tensor(xi(:,i+1),xi(:,l+1))'*Dxi(:,k)*Tensor(xi(:,i+1),xi(:,j+1),id2(:,r),xi(:,l+1));
                end
            end

            ABbasisS3Phi_Delta(:,k,j,r)=dummy/norm(dummy);

            
        end
    end
    r=2;
    for k=1:size(Dxi,2)
        for i=1:d-1
            dummy=zeros(2*d^3,1);
            for j=1:d-1
                for l=1:d-1
                    dummy=dummy+Tensor(xi(:,j+1),xi(:,l+1))'*Dxi(:,k)*Tensor(xi(:,i+1),xi(:,j+1),id2(:,r),xi(:,l+1));
                end
            end
            %    loccount=loccount+1;
            ABbasisS3Phi_Delta(:,k,i,r)=dummy/norm(dummy);
            %    R3(:,loccount)=ABbasisS3Phi_Delta(:,k,i,r);
        end
    end

    %% Copy of S3
    for k=1:size(Dxi,2)
        for i=1:d-1
            loccount=loccount+1;
            R3(:,loccount)=(ABbasisS3Phi_Delta(:,k,i,1)+ABbasisS3Phi_Delta(:,i,k,2))/norm(ABbasisS3Phi_Delta(:,k,i,1)+ABbasisS3Phi_Delta(:,i,k,2));
        end
    end



    %% Copy of S1
    for k=1:size(Dxi,2)
        for i=1:d-1
            loccount=loccount+1;
            R3(:,loccount)=(ABbasisS3Phi_Delta(:,k,i,1)-ABbasisS3Phi_Delta(:,i,k,2))/norm(ABbasisS3Phi_Delta(:,k,i,1)-ABbasisS3Phi_Delta(:,i,k,2));
        end
    end





    %S3Phi_Psi

    for k=1:size(Psiu,2)
        for j=2:d
            r=1;
            dummy=zeros(2*d^3,1);
            for i=2:d
                for l=2:d
                    dummy=dummy+Tensor(xi(:,i),xi(:,l))'*Psiu(:,k)*Tensor(xi(:,i),xi(:,j),id2(:,r),xi(:,l));
                end
            end
            loccount=loccount+1;
            ABbasisS3Phi_Delta(:,k,j,r)=dummy;
            R3(:,loccount)=ABbasisS3Phi_Delta(:,k,j,r);
        end
    end
    for k=1:size(Psiu,2)
        for i=2:d
            r=2;
            dummy=zeros(2*d^3,1);
            for j=2:d
                for l=2:d
                    dummy=dummy+Tensor(xi(:,j),xi(:,l))'*Psiu(:,k)*Tensor(xi(:,i),xi(:,j),id2(:,r),xi(:,l));
                end
            end
            loccount=loccount+1;
            ABbasisS3Phi_Delta(:,k,i,r)=dummy;
            R3(:,loccount)=ABbasisS3Phi_Delta(:,k,i,r);
        end
    end





    %% Phi Tensor Phi



    % Eplus
    for l=1:d-1
        for k=1:d-1
            loccount=loccount+1;
            ABbasisEplus(:,k,l)=(Tensor(XbasisPhi(:,k,1),YBbasisPhi(:,l,1))+Tensor(XbasisPhi(:,l,2),YBbasisPhi(:,k,2)))/norm(Tensor(XbasisPhi(:,k,1),YBbasisPhi(:,l,1))+Tensor(XbasisPhi(:,l,2),YBbasisPhi(:,k,2)));
            R3(:,loccount)=ABbasisEplus(:,k,l);
        end
    end



    % Eminus
    for l=1:d-1
        for k=1:d-1
            loccount=loccount+1;
            ABbasisEminus(:,k,l)=(Tensor(XbasisPhi(:,k,1),YBbasisPhi(:,l,1))-Tensor(XbasisPhi(:,l,2),YBbasisPhi(:,k,2)))/norm(Tensor(XbasisPhi(:,k,1),YBbasisPhi(:,l,1))-Tensor(XbasisPhi(:,l,2),YBbasisPhi(:,k,2)));
            R3(:,loccount)=ABbasisEminus(:,k,l);
        end
    end


    % FW
    loccount2=1;
    for k=1:size(Wbeta,2)
        ABbasisFWAbstract(:,loccount2)=Swap(Tensor(id2(:,1),id2(:,1),Wbeta(:,k)),[2 3],[2 2 d d])/norm(Tensor(id2(:,1),id2(:,1),Wbeta(:,k)));
        loccount2=loccount2+1;
    end
    for k=1:size(Wbeta,2)
        ABbasisFWAbstract(:,loccount2)=Swap(Tensor(id2(:,2),id2(:,2),Wbeta(:,k)),[2 3],[2 2 d d])/norm(Tensor(id2(:,2),id2(:,2),Wbeta(:,k)));
        loccount2=loccount2+1;
    end

    TransformOperator1=0;
    TransformOperator2=0;
    for i=2:d
        TransformOperator1=TransformOperator1+Tensor(xi(:,i),xi(:,1))*Tensor(id2(:,1),xi(:,i))';
        TransformOperator2=TransformOperator2+Tensor(xi(:,1),xi(:,i))*Tensor(id2(:,2),xi(:,i))';
    end
    TransformOperator=Tensor(TransformOperator1+TransformOperator2,id2,id); % in order to go from the (e_1\otimes\xi_i,e_2\otimes\xi_i) to (\xi_i\otimes\delta^+, \delta^+\otimes\xi_i)
    
    for k=1:size(ABbasisFWAbstract,2)
        ABbasisFW(:,k)=TransformOperator*ABbasisFWAbstract(:,k);
        loccount=loccount+1;
        R3(:,loccount)=ABbasisFW(:,k);
    end



    % FTplus
    ABbasisFTplusAbstract=Swap(Tensor(id2(:,1),id2(:,1),Tt)+Tensor(id2(:,2),id2(:,2),Tt),[2 3],[2 2 d d])/norm(Tensor(id2(:,1),id2(:,1),Tt)+Tensor(id2(:,2),id2(:,2),Tt));


    ABbasisFTplus=TransformOperator*ABbasisFTplusAbstract;
    loccount=loccount+1;
    R3(:,loccount)=ABbasisFTplus;

    % FTminus
    ABbasisFTminusAbstract=Swap(Tensor(id2(:,1),id2(:,1),Tt)-Tensor(id2(:,2),id2(:,2),Tt),[2 3],[2 2 d d])/norm(Tensor(id2(:,1),id2(:,1),Tt)-Tensor(id2(:,2),id2(:,2),Tt));

    ABbasisFTminus=TransformOperator*ABbasisFTminusAbstract;
    loccount=loccount+1;
    R3(:,loccount)=ABbasisFTminus;

    % FDelta
    loccount2=1;
    for l=2:-1:1
        for k=1:size(Dxi,2)

            ABbasisFDeltaAbstract(:,loccount2)=Swap(Tensor(id2(:,l),id2(:,l),Dxi(:,k)),[2 3],[2 2 d d])/norm(Tensor(id2(:,l),id2(:,l),Dxi(:,k)));
            loccount2=loccount2+1;
        end
    end


    for k=1:size(ABbasisFDeltaAbstract,2)
        ABbasisFDelta(:,k)=TransformOperator*ABbasisFDeltaAbstract(:,k);
        loccount=loccount+1;
        R3(:,loccount)=ABbasisFDelta(:,k);
    end


    % Fpsi
    ABbasisFpsiAbstract=[];
    loccount2=1;
    for l=2:-1:1
        for k=1:size(Psiu,2)
            ABbasisFpsiAbstract(:,loccount2)=Swap(Tensor(id2(:,l),id2(:,l),Psiu(:,k)),[2 3],[2 2 d d])/norm(Tensor(id2(:,l),id2(:,l),Psiu(:,k)));
            loccount2=loccount2+1;
        end
    end
    for k=1:size(ABbasisFpsiAbstract,2)
        ABbasisFpsi(:,k)=TransformOperator*ABbasisFpsiAbstract(:,k);
        loccount=loccount+1;
        R3(:,loccount)=ABbasisFpsi(:,k);
    end



    R=blkdiag(1,R1,R2,R3);




    p1=R(:,1);
    p2=R(:,2);
    p3=R(:,3:2+2*(d-1));
    p4=R(:,2+2*(d-1)+1:2+2*(d-1)+(d-1)^2);
    p5=R(:,3+2*(d-1)+(d-1)^2);
    p6=R(:,4+2*(d-1)+(d-1)^2);
    p7=R(:,5+2*(d-1)+(d-1)^2:4+4*(d-1)+(d-1)^2);
    p8=R(:,5+4*(d-1)+(d-1)^2);
    p9=R(:,6+4*(d-1)+(d-1)^2);
    p10=R(:,7+4*(d-1)+(d-1)^2:6+6*(d-1)+(d-1)^2); % t phi
    p11=R(:,7+6*(d-1)+(d-1)^2:6+8*(d-1)+(d-1)^2); %phi t
    p12=R(:,7+8*(d-1)+(d-1)^2:6+8*(d-1)+2*(d-1)^2);%phi S2
    p13=R(:,7+8*(d-1)+2*(d-1)^2:6+8*(d-1)+3*(d-1)^2);
    p14=R(:,7+8*(d-1)+3*(d-1)^2:6+10*(d-1)+3*(d-1)^2);

    %p15=R(:,7+10*(d-1)+3*(d-1)^2:6+10*(d-1)+3*(d-1)^2+2*(d-1)^3);%S3 phi
    p15W=R(:,7+10*(d-1)+3*(d-1)^2:6+10*(d-1)+3*(d-1)^2+(d-1)^2*(d-2));%S3 phi W
    p15t=R(:,7+10*(d-1)+3*(d-1)^2+(d-1)^2*(d-2):6+10*(d-1)+3*(d-1)^2+(d-1)^2*(d-2)+2*(d-1));%S3 phi
    p15d1=R(:,7+10*(d-1)+3*(d-1)^2+(d-1)^2*(d-2)+2*(d-1):6+10*(d-1)+3*(d-1)^2+(d-1)^2*(d-2)+2*(d-1)+(d-1)^2);%S3 delta 1
    p15d2=R(:,7+10*(d-1)+3*(d-1)^2+(d-1)^2*(d-2)+2*(d-1)+(d-1)^2:6+10*(d-1)+3*(d-1)^2+(d-1)^2*(d-2)+2*(d-1)+2*(d-1)^2);%S3 delta 2
    p15p=R(:,7+10*(d-1)+3*(d-1)^2+(d-1)^2*(d-2)+2*(d-1)+2*(d-1)^2:6+10*(d-1)+3*(d-1)^2+2*(d-1)^3);%S3 psi

    p16=R(:,7+10*(d-1)+3*(d-1)^2+2*(d-1)^3:6+10*(d-1)+4*(d-1)^2+2*(d-1)^3);
    p17=R(:,7+10*(d-1)+4*(d-1)^2+2*(d-1)^3:6+10*(d-1)+5*(d-1)^2+2*(d-1)^3);
    p18=R(:,7+10*(d-1)+5*(d-1)^2+2*(d-1)^3:6+10*(d-1)+5*(d-1)^2+2*(d-1)^3+2*(d-1)*(d-2)/2);%6+10*(d-1)+5*(d-1)^2+2*(d-1)^3+(d-1)*(d-2)/2);
                                                                                           %p19=R(:,7+10*(d-1)+5*(d-1)^2+2*(d-1)^3+(d-1)*(d-2)/2:);
    p20=R(:,7+10*(d-1)+5*(d-1)^2+2*(d-1)^3+2*(d-1)*(d-2)/2);
    p21=R(:,8+10*(d-1)+5*(d-1)^2+2*(d-1)^3+2*(d-1)*(d-2)/2);
    p22=R(:,9+10*(d-1)+5*(d-1)^2+2*(d-1)^3+2*(d-1)*(d-2)/2:8+12*(d-1)+5*(d-1)^2+2*(d-1)^3+2*(d-1)*(d-2)/2);%8+11*(d-1)+5*(d-1)^2+2*(d-1)^3+2*(d-1)*(d-2)/2);
                                                                                                           %p23=R(:,9+11*(d-1)+5*(d-1)^2+2*(d-1)^3+2*(d-1)*(d-2)/2:);
    p24=R(:,9+12*(d-1)+5*(d-1)^2+2*(d-1)^3+2*(d-1)*(d-2)/2:8+12*(d-1)+5*(d-1)^2+2*(d-1)^3+2*(d-1)*(d-2)/2+2*d*(d-3)/2);%8+12*(d-1)+5*(d-1)^2+2*(d-1)^3+2*(d-1)*(d-2)/2+d*(d-3)/2);
                                                                                                                       %p25=R(:,9+12*(d-1)+5*(d-1)^2+2*(d-1)^3+2*(d-1)*(d-2)/2+d*(d-3)/2:);



    NewR=[p1 p2 p5 p8 p20 p6 p9 p21 p3 p7 p10 p11 p15t p22 p14 p4 p12 p16 p15d1 p15d2 p13 p17 p18 p24 p15W p15p ];
    %p19 p22 p23 p24 p25 p15];

    function res = Tensor(varargin)
        res = [1];
        for iii = 1:length(varargin)
            res = kron(res, varargin{iii});
        end
    end

    function state = Swap(state, sys, dim)
        sys1 = length(dim) + 1 - sys(1);
        sys2 = length(dim) + 1 - sys(2);
        dim = fliplr(dim(:)');
        origdim = size(state);
        state = reshape(state, dim);
        sysperm = 1:length(dim);
        sysperm(sys2) = sys1;
        sysperm(sys1) = sys2;
        state = permute(state, sysperm);
        state = reshape(state, origdim);
    end

end
