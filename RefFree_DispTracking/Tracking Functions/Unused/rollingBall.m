% Rolling ball background subtraction
function filtImg = rollingBall(img,radius)
    % Structuring element, according to Sternberg 1983, should be a 3D
    % sphere whose diameter is considerably wider than any peaks yet small
    % enough to follow the smooth contours of the changing background
    % intensity. 

    % The method is unclear in the Sternberg paper, but this discussion
    % suggests that the Sternberg method is simply a top-hat transform with
    % a spherical structuring element. Further, there is no clear or
    % distinct advantage to using a spherical structuring element instead
    % of a simple disk for our purposes. In fact, a disk structring element
    % is arguably preferable due to faster performance. Therefore, in this
    % code, we will use a disk-shaped structuring element to perform the
    % "rolling ball" background subtraction. The code for performing the
    % rolling ball operation is left commented out in this code for
    % potential future use.
    %
    % Discussion:
    % http://dsp.stackexchange.com/questions/10597/uneven-background-subtraction-rolling-ball-vs-disk-tophat
    
    % The size of each dot in the gel is about 3-5 pixels. Therefore, the
    % specified radius should be greater than about 3.

% % % % %     % Rolling ball structuring element
% % % % %     se = strel('ball',radius,radius);

    se = strel('disk',radius);
    filtImg = imtophat(img,se);
%     imshow(filtImg,[])
end