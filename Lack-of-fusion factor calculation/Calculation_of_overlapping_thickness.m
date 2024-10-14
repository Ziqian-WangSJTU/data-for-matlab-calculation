clc;
clear all;
syms x y z tc R alpha0

%Set the index for calculated materials and number of points selected in the molten pool
number_of_material=[9 29];

%ambient temperatureambient temperature
T0=300;

%constant values for calculating incubation time
alpha1=0.45;
theta=40*pi/180;
ftheta=0.25*(2-3*cos(theta)+cos(theta)^3);
Sf=15.4388282915675; 
ar=10^-9;
da=2.55302157586516E-10;
D=0.000000002;

%laser power
P=[250 225 200];    

%scanning speed
V=[0.16 0.36 0.64];    

%physical properties for materials
material_property=readtable('F:\Printability prediction\matlab simulation\simulation required data\data_for_matlab.xlsx');

%molten pool size
opts=detectImportOptions('F:\Printability prediction\matlab simulation\simulation required data\lack of fusion simulation.xlsx');
opts.Sheet='molten pool size 2';
molten_pool=readtable('F:\Printability prediction\matlab simulation\simulation required data\lack of fusion simulation.xlsx',opts);

for powerindex=1:length(P)
    disp(P(powerindex))
    for velosityindex=1:length(V)
        disp(V(velosityindex))
        answer=[];
        for percentage_of_Ti=number_of_material
            molten_pool_size=table2array(molten_pool(:,['r_',num2str(percentage_of_Ti)]));
            hatching_space=table2array(molten_pool(:,['h_',num2str(percentage_of_Ti)]));
            powder_thichness=table2array(molten_pool(:,['d_',num2str(percentage_of_Ti)]));
            TL=table2array(material_property(percentage_of_Ti,'TL'));
            TS=table2array(material_property(percentage_of_Ti,'TS'));
            TC=table2array(material_property(percentage_of_Ti,'TC'));
            k=table2array(material_property(percentage_of_Ti,'k'));
            rou=table2array(material_property(percentage_of_Ti,'rou'));
            Cp=table2array(material_property(percentage_of_Ti,'Cp'));
            alpha=k/(rou*Cp);
            T=T0+alpha0.*P(powerindex).*(2*pi*k.*(x.^2+y.^2+z.^2).^0.5).^(-1).*exp(-V(velosityindex).*(x+(x.^2+y.^2+z.^2).^0.5)/(2.*alpha));   

            %adjusting the value of alpha0
            alpha0_test_value=0.25;
            temperature_in_x_z_plane=subs(T,{y,alpha0},{0,alpha0_test_value})-TL;
            temperary_gradient=gradient(temperature_in_x_z_plane,[x z]);
            dztodx=-temperary_gradient(1)/temperary_gradient(2);
            options=optimoptions('fsolve','Display','none');
            options.MaxFunctionEvaluations=600;
            Funtion_of_x_z=matlabFunction([temperature_in_x_z_plane;dztodx],'Vars',[x z]);
            Function_of_h=@(h) Funtion_of_x_z(h(1),h(2));
            molten_pool_radius=fsolve(Function_of_h,[0 -0.000001],options);
            alpha_left=0;
            alpha_right=1;
            effective_radius=molten_pool_size((powerindex-1)*3+velosityindex); 
            while abs(abs(molten_pool_radius(2))*10^6-effective_radius)>0.1 
                if abs(molten_pool_radius(2))*10^6>effective_radius
                    alpha_right=alpha0_test_value;
                else
                    alpha_left=alpha0_test_value;
                end
                alpha0_test_value=(alpha_left+alpha_right)/2;
                temperature_in_x_z_plane=subs(T,{y,alpha0},{0,alpha0_test_value})-TL;
                temperary_gradient=gradient(temperature_in_x_z_plane,[x z]);
                dztodx=-temperary_gradient(1)/temperary_gradient(2);
                Funtion_of_x_z=matlabFunction([temperature_in_x_z_plane;dztodx],'Vars',[x z]);
                Function_of_h=@(h) Funtion_of_x_z(h(1),h(2));
                molten_pool_radius=fsolve(Function_of_h,[0 -0.00001],options);
            end

            %Find the minimum point of the isothermal interface at TC
            T_function=subs(T,{y,alpha0},{0,alpha0_test_value})-TC;
            temperary_gradient=gradient(T_function,[x z]);
            T_dztodx=-temperary_gradient(1)/temperary_gradient(2);
            Funtion_of_x_z=matlabFunction([T_function;T_dztodx],'Vars',[x z]);
            Function_of_h=@(h) Funtion_of_x_z(h(1),h(2));
            minium_point_position=fsolve(Function_of_h,[0 -0.00001],options);
            
            line=subs(T,{x y alpha0},{minium_point_position(1) (hatching_space((powerindex-1)*3+velosityindex)*10^-6)/2 alpha0_test_value})-TC;
            answer=[answer double(vpasolve(line,[0 0.001]))*10^6];
        end         
        writematrix(answer,'Overlapping thickness.xlsx',"WriteMode","append");
    end
end
