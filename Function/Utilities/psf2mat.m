function [P] = psf2mat(p,sz,sf,shift)
% Aurthor: Ye Fei
% E-mail: feiye@njust.edu.cn
% Date: April 26, 2019  

    if nargin == 3
        shift=0;
    end
	if size(p,2)==1 && size(p,1)==1  && length(sz)>2
		error('PSF is not a matrix !!!');
    end
    if mod(size(p,1),2)~=0 && mod(size(p,2),2)~=0
        rotp=rot90(p,2);
        sz_p1=size(p,1);  sz_p2=size(p,2);
        for j=1:sz(2)
            if j<= sz_p2-ceil(sz_p2/2)+1
                if j==1
                    k=floor(sz_p2/2)+j;
                    r{j}=[rotp(ceil(sz_p1/2):end,k)' zeros(1,sz(1)-sz_p1) rotp(1:floor(sz_p1/2),k)'];
                    c{j}=[p(ceil(sz_p1/2):end,k)' zeros(1,sz(1)-sz_p1) p(1:floor(sz_p1/2),k)'];
                else
                    k=floor(sz_p2/2)+j;
                    r{j}=[rotp(ceil(sz_p1/2):end,k)' zeros(1,sz(1)-sz_p1) rotp(1:floor(sz_p1/2),k)'];
                    c{sz(2)+2-j}=[p(ceil(sz_p1/2):end,k)' zeros(1,sz(1)-sz_p1) p(1:floor(sz_p1/2),k)'];
                end
            elseif j>sz_p2-ceil(sz_p2/2)+1 && j<=sz(2)-floor(sz_p2/2)
                r{j}=zeros(1,sz(1));
                c{sz(2)+2-j}=zeros(1,sz(1));
            elseif j>sz(2)-floor(sz_p2/2)
                k=j-sz(2)+floor(sz_p2/2);
                r{j}=[rotp(ceil(sz_p1/2):end,k)' zeros(1,sz(1)-sz_p1) rotp(1:floor(sz_p1/2),k)'];
                c{sz(2)+2-j}=[p(ceil(sz_p1/2):end,k)' zeros(1,sz(1)-sz_p1) p(1:floor(sz_p1/2),k)'];
            end
        end
        for j=1:length(c)
            if j==1
%                 u{j}=sparse(toeplitz(c{j},r{j}));
                t=sparse(toeplitz(c{j},r{j}));
                u{j}=t(1+shift:sf:end,:);
%                 u{j}=t(1+shift:sf:end,1+shift:sf:end);
            else
%                 u{sz(2)+2-j}=sparse(toeplitz(c{j},r{j}));
                t=sparse(toeplitz(c{j},r{j}));
                u{sz(2)+2-j}=t(1+shift:sf:end,:);
%                 u{sz(2)+2-j}=t(1+shift:sf:end,1+shift:sf:end);
            end
        end
        s=speye(length(c));
%         s=t(1+shift:sf:end,:);
        P = sparse(zeros(prod(sz)/sf^2,prod(sz)));
        for j = 1:length(c)
            S = sparse(gallery('circul',s(j,:)')');
            S=S(1+shift:sf:end,:);
            P = P + kron(S,u{j});
        end
        
    else
        rotp=rot90(p,2);
        sz_p1=size(p,1);  sz_p2=size(p,2);
        for j=1:sz(2)
            if j<=sz_p2/2+1
                if j==1
                    k=sz_p2/2;
                    r{j}=[rotp(sz_p1/2:end,k)' zeros(1,sz(1)-sz_p1) rotp(1:sz_p1/2-1,k)'];
                    c{j}=[p(sz_p1/2+1:end,sz_p2+1-k)' zeros(1,sz(1)-sz_p1) p(1:sz_p1/2,sz_p2+1-k)'];
                else
                    k=sz_p2/2+j-1;
                    r{j}=[rotp(sz_p1/2:end,k)' zeros(1,sz(1)-sz_p1) rotp(1:sz_p1/2-1,k)'];
                    c{j}=[p(sz_p1/2+1:end,sz_p2+1-k)' zeros(1,sz(1)-sz_p1) p(1:sz_p1/2,sz_p2+1-k)'];
                end
            elseif j>sz_p2/2+1 && j<=sz(2)-sz_p2/2+1
                r{j}=zeros(1,sz(1));
                c{j}=zeros(1,sz(1));
            elseif j>sz(2)-sz_p2/2+1
                k=j-sz(2)+sz_p2/2-1;
                r{j}=[rotp(sz_p1/2:end,k)' zeros(1,sz(1)-sz_p1) rotp(1:sz_p1/2-1,k)'];
                c{j}=[p(sz_p1/2+1:end,sz_p2+1-k)' zeros(1,sz(1)-sz_p1) p(1:sz_p1/2,sz_p2+1-k)'];
            end
        end
        for j=1:length(c)
            if j==1
%                 u{j}=sparse(toeplitz(c{j},r{j}));
                t=sparse(toeplitz(c{j},r{j}));
                u{j}=t(1+shift:sf:end,:);
            else
%                 u{sz(2)+2-j}=sparse(toeplitz(c{j},r{j}));
                t=sparse(toeplitz(c{j},r{j}));
                u{sz(2)+2-j}=t(1+shift:sf:end,:);
            end
        end
        s=speye(length(c));
        P = sparse(zeros(prod(sz)/sf^2,prod(sz)));
        for j = 1 : length(c)
            S = sparse(gallery('circul',s(j,:)')');
            S=S(1+shift:sf:end,:);
            P = P + kron(S,u{j});
        end
    end
    P=P';
end