function display_inconsistency(PB,x,x_old)

if nargin==2
    x_old = [];
end

if isempty(x_old)
    q_old = [];
else
    q_old = get_q(PB,x_old);
end

q = get_q(PB,x);
[qsort,i] = sort(abs(q),'descend');
if length(i)>4
    qbig = abs(q) >= abs(q(i(4)));
else
    qbig = false(size(q));
end


for nq=1:PB.NQ
    inq1 = PB.L(1,nq);
    inq2 = PB.L(2,nq);
    v1 = PB.x2v(inq1,1);
    v2 = PB.x2v(inq2,1);
    if qbig(nq)
        sbig = '   !!!';
    else
        sbig = '';
    end
    disp(['Link #' num2str(nq) '; ' PB.varnames{v1} ' <--> ' PB.varnames{v2} sbig]);
    
    sold1 = '';
    if ~isempty(x_old)
        sold1 = [num2str(x_old(inq1)) ' --> '];
    end
    sold2 = '';
    if ~isempty(x_old)
        sold2 = [num2str(x_old(inq2)) ' --> '];
    end  
    disp(['    ' sold1 num2str(x(inq1)) '   //   ' sold2 num2str(x(inq2))]);
    
end