function seg_edge = ind2seg( ind )
% Input index "ind" marking segments of a vector. 0 for none segment.
% Return an m*2 matrix of index "seg_edge" specifying segment edges 
% where m is the number of segments.
seg_edge = reshape(find(diff([false;logical(ind(:));false])),2,[])';
seg_edge(:,2) = seg_edge(:,2)-1;
end