function [vec,val] = sortEig(V,order)
% wrapper for eig that sorts the eigenvectors and eigenvalues of V
% according to order. Vec and val are returned as tensors

[vec,val] = eig(V);
val = diag(val);
switch order
    case 'ascend'
        [val,idx] = sort(val,'ascend');
        tmp = vec;
        vec(:,1) = tmp(:,idx(1));
        vec(:,2) = tmp(:,idx(2));
        
        if numel(idx) > 2
            vec(:,3) = tmp(:,idx(3));
        end
        
    case 'descend'
        [val,idx] = sort(val,'descend');
        tmp = vec;
        vec(:,1) = tmp(:,idx(1));
        vec(:,2) = tmp(:,idx(2));
        if numel(idx) > 2
            vec(:,3) = tmp(:,idx(3));
        end        
    otherwise
        error('Please specify order as ascend or descend');
end
val = diag(val);

end