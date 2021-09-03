function flg = func_is_covered_by_circles(x_points, y_points, circles)
    tol = 1e-10;

    n_circle = size(circles,1);
    n_anchor = size(x_points, 1);
    
    x_circles = circles(:,1);
    y_circles = circles(:,2);
    r_circles = circles(:,3);
      
    mat_x_circle = repmat(x_circles', [n_anchor, 1]);
    mat_y_circle = repmat(y_circles', [n_anchor, 1]);
    mat_r_circle = repmat(r_circles', [n_anchor, 1]);
    
    mat_x_anchor = repmat(x_points, [1, n_circle]);
    mat_y_anchor = repmat(y_points, [1, n_circle]);
    
    mat_is_covered = (((mat_x_circle - mat_x_anchor).^2 + (mat_y_circle - mat_y_anchor).^2) <= (mat_r_circle.^2 - tol.*mat_r_circle));
    
    flg = any(mat_is_covered, 2);
end