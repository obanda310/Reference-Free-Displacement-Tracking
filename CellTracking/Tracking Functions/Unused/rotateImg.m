% Input image should be imregionalmax mask
function rotatedImg = rotateImg(img)
    [h,theta,rho] = hough(img);
    p = houghpeaks(h);
    lines = houghlines(img,theta,rho,p);
    % Create vector using the first line detected by houghlines by
    % subtracting the start-point coordinate from the end-point coordinate.
    % Note, z-dimension coordinate is zero.
    fp2 = [lines(1).point2,0];
    fp1 = [lines(1).point1,0];
    featureVector = fp2-fp1;
    % Create reference horizontal line vector that spans the entire width
    % of the image
    rp2 = [size(img,2),1,0];
    rp1 = [1,1,0];
    referenceVector = rp2-rp1;
    angle = atan2( ...
        norm(cross(featureVector,referenceVector)), ...
        dot(featureVector,referenceVector));
    % Output image rotated by calculated angle
    rotatedImg = imrotate(img,angle);
    
%     % Unecessary display of Hough-transformed image and detected lines    
%     imshow(imadjust(mat2gray(h)),'XData',theta,'YData',rho,...
%       'InitialMagnification','fit');
%     title('Hough transform of mask 10');
%     xlabel('\theta'), ylabel('\rho');
%     axis on, axis normal, hold on;
%     
%     figure;
%     imshow(img)
%     hold on
%     for i = 1:length(lines)
%         xy = [lines(i).point1; lines(i).point2];
%         plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green')
%     end
%     hold off
end