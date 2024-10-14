clc;
clear all;
syms x y z tc R alpha0

%Set the index for calculated materials and number of points selected in the molten pool
number_of_material=[9 29];
number_of_point_y=21;
number_of_point_z=19;

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

%constant value in KGT dendritic growth model
A=1042;
n=3.1;

%physical properties for materials
material_property=readtable('F:\Printability prediction\matlab simulation\simulation required data\data_for_matlab.xlsx');

%actual molten pool size measured
molten_pool_size=readtable('F:\Printability prediction\matlab simulation\simulation required data\molten pool size.xlsx');

for powerindex=1:1:length(P)
    disp(P(powerindex))
    for velosityindex=1:1:length(V)
        disp(V(velosityindex))
        for percentage_of_Ti=number_of_material
            filename=[num2str(((percentage_of_Ti+1)/20)),' ',num2str(P(powerindex)),' ',num2str(V(velosityindex)),'.xlsx'];  %saving path for output files
            TL=table2array(material_property(percentage_of_Ti,'TL'));
            TS=table2array(material_property(percentage_of_Ti,'TS'));
            TC=table2array(material_property(percentage_of_Ti,'TC'));
            k=table2array(material_property(percentage_of_Ti,'k'));
            rou=table2array(material_property(percentage_of_Ti,'rou'));
            Cp=table2array(material_property(percentage_of_Ti,'Cp'));
            alpha=k/(rou*Cp);
            T=T0+alpha0.*P(powerindex).*(2*pi*k.*(x.^2+y.^2+z.^2).^0.5).^(-1).*exp(-V(velosityindex).*(x+(x.^2+y.^2+z.^2).^0.5)/(2.*alpha));   
            y0=table2array(material_property(percentage_of_Ti,'y0'));
            x0=table2array(material_property(percentage_of_Ti,'x0'));
            A1=table2array(material_property(percentage_of_Ti,'A1'));
            t1=table2array(material_property(percentage_of_Ti,'t1'));

            %adjusting the value of alpha0
            alpha0_test_value=0.25;
            temperature_in_x_z_plane=subs(T,{y,alpha0},{0,alpha0_test_value})-TL;%
            temperary_gradient=gradient(temperature_in_x_z_plane,[x z]);
            dztodx=-temperary_gradient(1)/temperary_gradient(2); %find the minimum point of molten pool
            options=optimoptions('fsolve','Display','none');
            options.MaxFunctionEvaluations=600;
            Funtion_of_x_z=matlabFunction([temperature_in_x_z_plane;dztodx],'Vars',[x z]);
            Function_of_h=@(h) Funtion_of_x_z(h(1),h(2));
            molten_pool_radius=fsolve(Function_of_h,[0 -0.000001],options);
            alpha_left=0;
            alpha_right=1;
            effective_radius=table2array(molten_pool_size((powerindex-1)*3+velosityindex,(percentage_of_Ti+11)/20));
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

            %final temperature field in the molten pool
            terminal_temperature_field=subs(T,alpha0,alpha0_test_value);
            temperature_gradient_field=gradient(terminal_temperature_field,[x y z]);
            sheetname=['Ti ',num2str((percentage_of_Ti+1)/20)];
            
            sum_F=0;
            number_of_F=0;
            disp(percentage_of_Ti)

            %calculation for F
            for pointindex_in_z=1:9/(number_of_point_z-1):10
                value_of_z=-pointindex_in_z*10^-5;
                if abs(value_of_z)>abs(molten_pool_radius(2))
                    break
                end
                coolingrate_tmp=[];
                G_tmp=[];
                V_solidification_tmp=[];
                Al3Ti_tmp=[];
                phi_tmp=[];
                F_tmp=[];
                for pointindex_in_y=0:10/(number_of_point_y-1):10
                    value_of_y=pointindex_in_y*10^-5;
                    if (value_of_y^2+value_of_z^2)^0.5>abs(molten_pool_radius(2))
                        break
                    end
                    T_x_curve=subs(terminal_temperature_field,{x y z},{-V(velosityindex)*tc value_of_y value_of_z});
                    position_of_T_max=vpasolve(diff(T_x_curve,tc),[0 Inf]);
                    tc_of_TL_isothermal_face=double(vpasolve(T_x_curve-TL,[position_of_T_max Inf]));
                    tc_of_TC_isothermal_face=double(vpasolve(T_x_curve-TC,[position_of_T_max Inf]));
                    
                    coolingrate=abs(double(integral(matlabFunction(diff(T_x_curve,tc)),tc_of_TL_isothermal_face,tc_of_TC_isothermal_face)/(tc_of_TC_isothermal_face-tc_of_TL_isothermal_face)));
                    G_vector=subs(temperature_gradient_field,{x y z},{-tc_of_TC_isothermal_face*V(velosityindex) value_of_y value_of_z});
                    G=double(norm(G_vector));
                    V_solidification=V(velosityindex)*G_vector(1)/G;

                    Tr=(TL-coolingrate*tc)/TL;
                    tresidual=(TL-TC)/coolingrate;

                    %fited function for Xeff
                    xLeff=y0+A1*exp(((TL-coolingrate*tc)-x0)/t1);
                     

                    %calculation of τ
                    tao=7.2*8.314*ftheta/(1-cos(theta))*ar^4/(da^2.*xLeff).*Tr.*(D*Sf*(1-Tr).^2).^-1;
                    th=tc-pi.^2/6.*tao;
                    Al3Ti_Number_function=exp(-16*pi*alpha1^3*ftheta/(3*8.314)*Sf*(1-Tr)^(-2)/Tr)*th;
                    if subs(th,tc,tresidual)>0
                        tl=double(vpasolve(th,tc,[0 tresidual]));
                        average_number_of_Al3Ti=integral(matlabFunction(Al3Ti_Number_function),tl,tresidual)*10^25/(tresidual-tl);
                    else
                        if subs(diff(th,tc),tc,tresidual)>0
                            average_number_of_Al3Ti=10^9;
                        else
                            top_time=vpasolve(diff(th,tc),[0.000001 tresidual]);
                            if subs(th,tc,top_time)>0
                                tl=double(vpasolve(th,[0 top_time]));
                                th=double(vpasolve(th,[top_time tresidual]));
                                average_number_of_Al3Ti=integral(matlabFunction(Al3Ti_Number_function),tl,th)*10^25/(th-tl);
                            else
                                average_number_of_Al3Ti=10^9;
                            end
                        end
                    end

                    %calculation for φ
                    phi=1-exp(-4*pi*average_number_of_Al3Ti*(V_solidification*A)^(3/n)/(3*((n+1)*G)^3));

                    %calculation for F
                    G_full_collumnar=(-4*pi*average_number_of_Al3Ti/(3*log(1-0.0066)))^(1/3)/(n+1)*(A*R)^(1/n);
                    G_full_equaixed=(-4*pi*average_number_of_Al3Ti/(3*log(1-0.66)))^(1/3)/(n+1)*(A*R)^(1/n);
                    M=G*V_solidification;
                    assistant_line=M/R;
                    theta1=atan(M/(vpasolve(G_full_equaixed-assistant_line,R,[-Inf Inf]))^2);
                    theta2=atan(M/(vpasolve(G_full_collumnar-assistant_line,R,[-Inf Inf]))^2);
                    theta3=atan(G/V_solidification);
                    alpha_theta=theta3-theta1;
                    beta_theta=theta2-theta3;
                    if alpha_theta/beta_theta>-1
                        F=1.5/(alpha_theta/beta_theta+3)+0.25;
                    else
                        F=1.5/(alpha_theta/beta_theta-5)+0.25;
                    end

                    sum_F=sum_F+F;
                    number_of_F=number_of_F+1;

                    %record the values that need to be outputed
                    coolingrate_tmp=[coolingrate_tmp coolingrate];
                    G_tmp=[G_tmp G];
                    V_solidification_tmp=[V_solidification_tmp V_solidification];
                    Al3Ti_tmp=[Al3Ti_tmp average_number_of_Al3Ti];
                    phi_tmp=[phi_tmp phi];
                    F_tmp=[F_tmp F];
                end
                writematrix(double(coolingrate_tmp),filename,'Sheet',['cooling rate ',sheetname],'WriteMode','append')
                writematrix(double(G_tmp),filename,'Sheet',['temperature gradient ',sheetname],'WriteMode','append')
                writematrix(double(V_solidification_tmp),filename,'Sheet',['solidification speed ',sheetname],'WriteMode','append')
                writematrix(double(Al3Ti_tmp),filename,'Sheet',['Al3Ti number density ',sheetname],'WriteMode','append')
                writematrix(double(phi_tmp),filename,'Sheet',['phi ',sheetname],'WriteMode','append')
                writematrix(double(F_tmp),filename,'Sheet',['F ',sheetname],'WriteMode','append')
            end
            writematrix([double(sum_F/number_of_F) alpha0_test_value],filename,'Sheet',['F ',sheetname],'WriteMode','append')
        end         
    end
end
