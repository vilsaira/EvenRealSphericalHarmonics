classdef SH_animator < handle
    properties
        lmax = [];
        Ncoeffs = [];
        Npoints = [];
        DMat = [];
        coeffs = [];
        signal = [];
        x = [];
        y = [];
        z = [];
        sx = [];
        sy = [];
        sz = [];
        c = [];
        phi = [];
        theta = [];
        l_list = [];
        m_list = [];
    end
    methods
        function obj = SH_animator(lmax, Npoints)
            obj.lmax = lmax;
            obj.Npoints = Npoints;
            obj.Ncoeffs = (lmax+2)*(lmax+1)/2;
            obj.coeffs = zeros(obj.Ncoeffs,1);
            [obj.sx,obj.sy,obj.sz] = sphere(obj.Npoints);
            obj.calcDMat();
            [obj.signal, obj.c, obj.phi, obj.theta, obj.x, obj.y, obj.z] = ...
                obj.calcSignal(obj, obj.coeffs);
        end
        function obj = calcDMat(obj)
            [az, el, ~] = cart2sph(obj.sx,obj.sy,obj.sz);
            az = az(:); % az [-pi, pi]
            el = el(:); % el [-pi/2, pi/2]

            obj.theta = el+pi/2;    % colatitude angle [0, pi]
            obj.phi = az+pi;        % azimuthal angle [0, 2*pi]

            Vpoints = length(obj.phi);
            obj.DMat = zeros(Vpoints, obj.Ncoeffs);

            col = 0;
            for l = 0:2:obj.lmax % only even degrees
                Plm = legendre(l, cos(obj.theta)); % lenegdre calculates all m>=0 for a given l

                for m = -l:1:+l
                    col = col + 1;
                    obj.l_list(col) = l;
                    obj.m_list(col) = m;
                    absm = abs(m);
                    if l > 0
                        plm = reshape(Plm(absm+1,:,:), size(obj.phi));
                    else
                        plm = reshape(Plm(:,:), size(obj.phi));
                    end

                    a = (2*l+1)*factorial(l-absm);
                    b = 4*pi*factorial(l+absm);
                    C = sqrt(a/b);
                    Ylm = C.*plm.*exp(1i * absm * obj.phi);

                    SH_coeffs = real(Ylm);
                    if m < 0
                        SH_coeffs = imag(Ylm);
                    end
                    if m ~= 0
                        SH_coeffs = sqrt(2)*SH_coeffs;
                    end

                    obj.DMat(:, col) = SH_coeffs;
                end
            end
        end
    end
    methods (Static)
        function [signal, c, phi, theta, x, y, z] = calcSignal(obj, coeffs)
            signal = obj.DMat*coeffs;
            c = reshape(signal, size(obj.sx));
            phi = reshape(obj.phi, size(obj.sx));
            theta = reshape(obj.theta, size(obj.sx));
            [x, y, z] = sph2cart(phi, pi/2-theta, c);
        end
        function fh = draw(obj, inds, vals, type)
            coeffs = obj.coeffs;
            coeffs(inds) = vals;
            [obj.signal, obj.c, obj.phi, obj.theta, obj.x, obj.y, obj.z] = ...
                obj.calcSignal(obj, coeffs);
            fh = figure('Position', [0 0 1920 1080]);
            switch type
                case 1
                    s = surf(obj.x, obj.y, obj.z, obj.c,'edgecolor',[3, 160, 98]./255,...
                        'FaceLighting', 'gouraud', 'FaceColor', 'none');
                        lightangle(0,0);
                case 2
                    s = surf(obj.x, obj.y, obj.z, obj.c,'edgecolor','none',...
                        'FaceLighting', 'gouraud');
                        lightangle(0,0);
            end
            shading interp;
            lighting gouraud;
            axis equal;
            material dull;
            axis off;
            lh = light;
            set(gcf,'color','k');
            
            s.FaceLighting = 'gouraud';
            s.AmbientStrength = 0.4;
            s.DiffuseStrength = 0.5;
            s.SpecularStrength = 0.1;
            s.SpecularExponent = 10;
            s.BackFaceLighting = 'unlit';
            axis tight;
            InSet = get(gca, 'TightInset');
            set(gca, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])
            
            colormap hsv;
            camzoom(1.0);
            view(-15, -55);
            %% Animation
            dtheta = 3;
            lighttheta = 0;
            R = linspace(0,1,60);
            pd = makedist('HalfNormal','mu',0,'sigma',0.3);
            T = pdf(pd, R);
            T = T - min(T);
            T = T./max(T);
            T = 1-T;
%             plot(R,T);
          
            vidfile = VideoWriter('testmovie2.mp4', 'MPEG-4');
            vidfile.FrameRate = 60;
            open(vidfile);
            col = 0;
            frameN = 0;
            for l = 0:2:obj.lmax
                for m = -l:1:+l
                    col = col + 1;
                    if col < obj.Ncoeffs
                        coeffs = zeros(size(coeffs));
                        for ri = 1:length(R)
                            frameN = frameN + 1;
                            coeffs(col+1) = T(ri);
                            if frameN == 1
                                coeffs(col+1) = T(ri+1);
                            end
                            coeffs(col) = 1-T(ri);
                            [obj.signal, obj.c, obj.phi, obj.theta, obj.x, obj.y, obj.z] = ...
                                obj.calcSignal(obj, coeffs);
                            s.XData = obj.x;
                            s.YData = obj.y;
                            s.ZData = obj.z;
                            s.CData = obj.c;
                            lighttheta = lighttheta + dtheta;
                            camorbit(1,0,'camera');
%                             lightangle(lh,0,lighttheta);
                            F = getframe(fh);
                            writeVideo(vidfile, F);
                            pause(0.001)
                        end
                    end
                    disp(col / obj.Ncoeffs);
                end
            end
            close(vidfile)
        end
    end
end
