function contact_point = find_contact_point(leg_contour , landscape_table , r)
% When the leg_contour and the gound is overlap, return the contact point
% and the distance needed as the displacement of thew leg

    error_dis = 0.01;


    contact_area_index_1 = leg_contour.leg_1.contour.y <= interp1(landscape_table(1,:),landscape_table(2,:),leg_contour.leg_1.contour.x);
    contact_area_index_2 = leg_contour.leg_2.contour.y <= interp1(landscape_table(1,:),landscape_table(2,:),leg_contour.leg_2.contour.x);
    
    % find all the points of contour lower than the landscape
    
    leg_1_touch_gound = ~isempty( find(contact_area_index_1, 1) ) ;
    leg_2_touch_gound = ~isempty( find(contact_area_index_2, 1) ) ;


    if leg_1_touch_gound
        % define searching area of x of the landscape
        contact_area_x_1 = leg_contour.leg_1.contour.x( 1, contact_area_index_1 );

        contact_area_y_1 = interp1(landscape_table(1,:),landscape_table(2,:),contact_area_x_1);
% 
%         % find the lowest point
%         [min_dif_1,min_dif_index_1] = ...
%          min( leg_contour.leg_1.contour.y(1, contact_area_index_1) - interp1(landscape_table(1,:),landscape_table(2,:),contact_area_x_1));
% 
%         contact_point_x = contact_area_x_1(min_dif_index_1);
% 
%         contact_point_1 = [contact_point_x, interp1(landscape_table(1,:),landscape_table(2,:),contact_point_x)];
% 
% 
%         revise_vector_1 = (leg_contour.leg_1.center - contact_point_1) / norm(leg_contour.leg_1.center - contact_point_1);
%         revise_vector_1 = revise_vector_1 * min_dif_1;
% 
%         contact_point.point_1 = [contact_point_1 , revise_vector_1];
         

     %%%%%%%% solution 2
        landscape_points_1 = [contact_area_x_1 ; contact_area_y_1];   
        ditance_matrix_1 = pdist2 (leg_contour.leg_1.center , landscape_points_1'  );
        % ditance_matrix(i,j) 
        % means the distance from ith point in X to jth point in Y
        
        % Find the min_dis between the half circle center and the ground
        min_distance_1 = min(ditance_matrix_1(:));
        [i,j] = find(ditance_matrix_1 == min_distance_1 , 1);
        
        revise_vector_1 = (leg_contour.leg_1.center - landscape_points_1(:,j)')/norm(leg_contour.leg_1.center - landscape_points_1(:,j)');
        revise_vector_1 = revise_vector_1*(r - min_distance_1);
        
        contact_point.point_1.point = landscape_points_1(:,j)';
        contact_point.point_1.revise = revise_vector_1;

        %ditance between contact point and tip & toe
        contact_point.point_1.istip = norm(landscape_points_1(:,j)' - leg_contour.leg_1.tip ) < error_dis;
        contact_point.point_1.istoe = norm(landscape_points_1(:,j)' - leg_contour.leg_1.toe ) < error_dis;
 
    else
        contact_point.point_1.point = [];
        contact_point.point_1.revise = [];
    end

    if leg_2_touch_gound

        contact_area_x_2 = leg_contour.leg_2.contour.x( 1, contact_area_index_2 );
        contact_area_y_2 = interp1(landscape_table(1,:),landscape_table(2,:),contact_area_x_2);

%         [min_dif_2,min_dif_index_2] = ...
%          min( leg_contour.leg_2.contour.y(1, contact_area_index_2) - interp1(landscape_table(1,:),landscape_table(2,:),contact_area_x_2));
%         contact_point_x = contact_area_x_2(min_dif_index_2);
%         
%         contact_point_2 = [contact_point_x, interp1(landscape_table(1,:),landscape_table(2,:),contact_point_x) ];
%         
%         revise_vector_2 = (leg_contour.leg_2.center - contact_point_2) / norm(leg_contour.leg_2.center - contact_point_2);
%         revise_vector_2 = revise_vector_2 * min_dif_2;
%         
%         contact_point.point_2 = [contact_point_2 , revise_vector_2];


    %%%%%%%% solution 2
        landscape_points_2 = [contact_area_x_2 ; contact_area_y_2];   
        ditance_matrix_2 = pdist2 ( leg_contour.leg_2.center , landscape_points_2'  );
        % ditance_matrix(i,j) 
        % means the distance from ith point in X to jth point in Y
        
        min_distance_2 = min(ditance_matrix_2(:));
        [u,v] = find(ditance_matrix_2 == min_distance_2 , 1);
        revise_vector_2 = (leg_contour.leg_2.center - landscape_points_2(:,v)') / norm(leg_contour.leg_2.center - landscape_points_2(:,v)');
        revise_vector_2 = revise_vector_2 * (r - min_distance_2);
        contact_point.point_2.point = landscape_points_2(:,v)';
        contact_point.point_2.revise = revise_vector_2;

        %ditance between contact point and tip & toe
        contact_point.point_2.istip = norm(landscape_points_2(:,v)' - leg_contour.leg_2.tip ) < error_dis;
        contact_point.point_2.istoe = norm(landscape_points_2(:,v)' - leg_contour.leg_2.toe ) < error_dis;
    else
        contact_point.point_2.point = [];
        contact_point.point_2.revise = [];
    end


end



