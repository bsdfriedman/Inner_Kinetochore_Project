function create_IKK_pericen(chrom_rad_nm,arc_length_nm,num_DNA_beads)

% input check
if mod(num_DNA_beads,2) == 0
    error('num_DNA_beads must be an odd integer')
else
end

% Parameters
pin_2_cen = 1.5*10^-7;
cen_2_mt = 5*10^-8;
mt_radius = 1.25*10^-7;

% calculate theta to create the cen bead and the mt bead positions
theta = (2*pi)/16;

% preallocate
cen_coords = zeros([16 3]);

% place the cen coords
for z = 1:16
    cen_coords(z,1:3) = [1 0 0; 0 cos(theta*(z-1)) -sin(theta*(z-1)); 0 sin(theta*(z-1)) cos(theta*(z-1))]*[0; mt_radius; 0];
end

% create the mt coords
mt_coords = cen_coords;
mt_coords(:,1) = cen_2_mt;

% create the base vectors from which the pinned chromosome vectors will be based
pin_base_vecs = cell(16,1);
for z = 1:size(cen_coords,1)
    pin_base_vecs{z,1} = cat(2,-pin_2_cen,([cen_coords(z,2) cen_coords(z,3)] / norm(cen_coords(z,:)))*chrom_rad_nm*10^-9);
end

% calculate the rotation angle for the pinned DNA
theta = (arc_length_nm*10^-9/2)/(chrom_rad_nm*10^-9);

% use vector rotation and the base coords to calculate the pinned coords
pin_coords = zeros([size(pin_base_vecs,1)*2 3]);
for z = 1:size(pin_base_vecs,1)
    pos_new_1 = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)]*pin_base_vecs{z,1}';
    pos_new_2 = [1 0 0; 0 cos(-theta) -sin(-theta); 0 sin(-theta) cos(-theta)]*pin_base_vecs{z,1}';
    pin_coords(z*2-1,1:3) = pos_new_1(1:3);
    pin_coords(z*2,1:3) = pos_new_2(1:3);
end

% use vector rotation to make reference coords for centromeres to build the chains
cen_ref_coords = zeros([size(cen_coords,1)*2 3]);
for z = 1:size(cen_coords,1)
    pos_new_1 = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)]*cen_coords(z,1:3)';
    pos_new_2 = [1 0 0; 0 cos(-theta) -sin(-theta); 0 sin(-theta) cos(-theta)]*cen_coords(z,1:3)';
    cen_ref_coords(z*2-1,1:3) = pos_new_1(1:3);
    cen_ref_coords(z*2,1:3) = pos_new_2(1:3);
end

% preallocate the chain coords matrix
chain_coords = zeros([size(pin_base_vecs,1)*(num_DNA_beads-3) 3]);

% place the additional DNA beads between the pinned end and the centromeres
half_chain_beads = (num_DNA_beads-1)/2;
if norm(pin_coords(1,:)-cen_ref_coords(1,:)) >= half_chain_beads*10^-8
    % if the pin-cen distance is long enough for a straight line
    % create the line vector and evenly space the beads
    for z = 1:size(pin_base_vecs,1)
        in_vec = (cen_ref_coords(z*2-1,:)-pin_coords(z*2-1,:))/half_chain_beads;
        out_vec = (pin_coords(z*2,:)-cen_ref_coords(z*2,:))/half_chain_beads;
        % loop through and use vector addition to place the remaining beads
        for h = 1:half_chain_beads-1
            chain_coords(((z-1)*(half_chain_beads-1)*2)+h,1:3) = pin_coords(z*2-1,:) + in_vec*h;
            chain_coords(((z-1)*(half_chain_beads-1)*2)+(half_chain_beads-1)+h,1:3) = cen_ref_coords(z*2,:) + out_vec*h;
        end
    end
else
    % if the pin-cen distance is close enough to require a triangle
    % calculate the hypotenuse length and the height of the triangle
    hypot = half_chain_beads*10^-8/2;
    tri_height = sqrt(hypot^2 - ((norm(pin_coords(1,:)-cen_ref_coords(1,:)))/2)^2);
    for z = 1:size(pin_base_vecs,1)
        % create the unit vectors and find the vertex of the triangle
        unit_vec_1 = [0 pin_coords(z*2-1,2) pin_coords(z*2-1,3)] / norm([0 pin_coords(z*2-1,2) pin_coords(z*2-1,3)]);
        unit_vec_2 = [0 pin_coords(z*2,2) pin_coords(z*2,3)] / norm([0 pin_coords(z*2,2) pin_coords(z*2,3)]);
        tri_vert_1 = ((pin_coords(z*2-1,:)+cen_ref_coords(z*2-1,:))/2) + (unit_vec_1*tri_height);
        tri_vert_2 = ((pin_coords(z*2,:)+cen_ref_coords(z*2,:))/2) + (unit_vec_2*tri_height);
        
        % create vectors to plot the remaining beads
        if mod(half_chain_beads-1,2) == 1
            % if there is an odd number of beads
            % create the placement vectors
            up_vec_1 = (tri_vert_1 - pin_coords(z*2-1,:))/((half_chain_beads)/2);
            up_vec_2 = (tri_vert_2 - cen_ref_coords(z*2,:))/((half_chain_beads)/2);
            down_vec_1 = (cen_ref_coords(z*2-1,:) - tri_vert_1)/((half_chain_beads)/2);
            down_vec_2 = (pin_coords(z*2,:) - tri_vert_2)/((half_chain_beads)/2);
            
            % create the new beads
            for h = 1:(half_chain_beads)/2
                chain_coords(((z-1)*(half_chain_beads-1)*2)+h,1:3) = pin_coords(z*2-1,:) + up_vec_1*h;
                chain_coords(((z-1)*(half_chain_beads-1)*2)+(half_chain_beads-1)+h,1:3) = cen_ref_coords(z*2,:) + up_vec_2*h;
            end
            for h = 1+((half_chain_beads)/2):half_chain_beads-1
                chain_coords(((z-1)*(half_chain_beads-1)*2)+h,1:3) = tri_vert_1 + down_vec_1*(h-((half_chain_beads)/2));
                chain_coords(((z-1)*(half_chain_beads-1)*2)+(half_chain_beads-1)+h,1:3) = tri_vert_2 + down_vec_2*(h-((half_chain_beads)/2));
            end
        else
            % if there is an even number of beads
            % create the placement vectors
            up_vec_1 = (tri_vert_1 - pin_coords(z*2-1,:))/((half_chain_beads-1)/2);
            up_vec_2 = (tri_vert_2 - cen_ref_coords(z*2,:))/((half_chain_beads-1)/2);
            down_vec_1 = (cen_ref_coords(z*2-1,:) - tri_vert_1)/((half_chain_beads-1)/2);
            down_vec_2 = (pin_coords(z*2,:) - tri_vert_2)/((half_chain_beads-1)/2);
            
            % create the new beads
            for h = 1:(half_chain_beads-1)/2
                chain_coords(((z-1)*(half_chain_beads-1)*2)+h,1:3) = pin_coords(z*2-1,:) + up_vec_1*(h-0.5);
                chain_coords(((z-1)*(half_chain_beads-1)*2)+(half_chain_beads-1)+h,1:3) = cen_ref_coords(z*2,:) + up_vec_2*(h-0.5);
            end
            for h = 1+((half_chain_beads-1)/2):half_chain_beads-1
                chain_coords(((z-1)*(half_chain_beads-1)*2)+h,1:3) = tri_vert_1 + down_vec_1*(h-0.5-((half_chain_beads-1)/2));
                chain_coords(((z-1)*(half_chain_beads-1)*2)+(half_chain_beads-1)+h,1:3) = tri_vert_2 + down_vec_2*(h-0.5-((half_chain_beads-1)/2));
            end
        end
    end
end

% assign chromoShake parameters
mass_mass = 3.38889e-020; % mass of a standard bead
mass_multiplier = 10^10; % pinning mass
spring_rest = 1e-008; % spring distance at rest
spring_const = 0.226195; % standard spring constant
hinge_const = 4.0715e-012; % standard hinge constant

% intro for chromoShake
intro = sprintf('meta temperature_Celsius 25\r\nmeta viscosity_centiPoise 1\r\nmeta effective_damping_radius 8e-009\r\nmeta dna_modulus_gigaPascal 2\r\nmeta dna_radius_nanometers 0.6\r\nmeta damping_radius_factor 0.8\r\nstructure {\r\n  random_force 2.78554e-011\r\n  mass_damping 4.44973e+009\r\n  mass_radius 4.5e-009\r\n  time_step 2e-009\r\n  collision_spring_constant 0.0565487\r\n  spring_damping_factor 0\r\n  random_number_seed 42\r\n  color 1');

% open up the file
outfile_name = sprintf('IKK_pericentromere_r%d_a%d_b%d_s42.cfg',chrom_rad_nm,arc_length_nm,num_DNA_beads);
fid_out = fopen(outfile_name,'w');

% print out the intro
fprintf(fid_out,'%s\r\n',intro);

% print out the masses, springs, and hinges
for z = 1:size(pin_base_vecs,1)
    % print the first pinned end
    fprintf(fid_out,'  mass %d\t%.6g\t%.6g %.6g %.6g %d\r\n',(num_DNA_beads+1)*(z-1),mass_mass*mass_multiplier,pin_coords(z*2-1,1),pin_coords(z*2-1,2),pin_coords(z*2-1,3),2);
    
    % print the first arm
    for h = 1:half_chain_beads-1
        fprintf(fid_out,'  mass %d\t%.6g\t%.6g %.6g %.6g %d\r\n',((num_DNA_beads+1)*(z-1))+h,mass_mass,chain_coords((((half_chain_beads-1)*2)*(z-1))+h,1),chain_coords((((half_chain_beads-1)*2)*(z-1))+h,2),chain_coords((((half_chain_beads-1)*2)*(z-1))+h,3),1);
    end
    
    % print the centromere bead
    fprintf(fid_out,'  mass %d\t%.6g\t%.6g %.6g %.6g %d\r\n',(num_DNA_beads+1)*(z-1)+half_chain_beads,mass_mass,cen_coords(z,1),cen_coords(z,2),cen_coords(z,3),2);
    
    % print the second arm
    for h = 1:half_chain_beads-1
        fprintf(fid_out,'  mass %d\t%.6g\t%.6g %.6g %.6g %d\r\n',((num_DNA_beads+1)*(z-1))+half_chain_beads+h,mass_mass,chain_coords((((half_chain_beads-1)*2)*(z-1))+half_chain_beads-1+h,1),chain_coords((((half_chain_beads-1)*2)*(z-1))+half_chain_beads-1+h,2),chain_coords((((half_chain_beads-1)*2)*(z-1))+half_chain_beads-1+h,3),1);
    end
    
    % print the second pinned end
    fprintf(fid_out,'  mass %d\t%.6g\t%.6g %.6g %.6g %d\r\n',(num_DNA_beads+1)*(z)-2,mass_mass*mass_multiplier,pin_coords(z*2,1),pin_coords(z*2,2),pin_coords(z*2,3),2);
    
    % print the microtubule bead
    fprintf(fid_out,'  mass %d\t%.6g\t%.6g %.6g %.6g %d\r\n',(num_DNA_beads+1)*(z)-1,mass_mass*mass_multiplier,mt_coords(z,1),mt_coords(z,2),mt_coords(z,3),3);
    
    % create the DNA springs
    for h = ((num_DNA_beads+1)*(z-1)):(num_DNA_beads+1)*(z)-3
        fprintf(fid_out,'  spring %d %d %.1g %.6g\r\n',h,h+1,spring_rest,spring_const);
    end
    
    % create the DNA hinges
    for h = ((num_DNA_beads+1)*(z-1)):(num_DNA_beads+1)*(z)-4
        fprintf(fid_out,'  hinge %d %d %d %.6g\r\n',h,h+1,h+2,hinge_const);
    end
    
    % create the microtubule spring
    fprintf(fid_out,'  spring %d %d %.1g %.6g\r\n',(num_DNA_beads+1)*(z-1)+half_chain_beads,(num_DNA_beads+1)*(z)-1,cen_2_mt,spring_const);
    
end

% print out the end of the structure
fprintf(fid_out,'}\r\n');

fclose('all');

% EXTRA CODE FOR PLOTTING POINTS:

% scatter3(mt_coords(:,1),mt_coords(:,2),mt_coords(:,3));
% hold on;
% scatter3(cen_coords(:,1),cen_coords(:,2),cen_coords(:,3));
% scatter3(pin_coords(:,1),pin_coords(:,2),pin_coords(:,3));
% scatter3(chain_coords(:,1),chain_coords(:,2),chain_coords(:,3));
% hold off;
