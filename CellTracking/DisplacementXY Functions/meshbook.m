%% Code used to work, but needs to be updated. Code was migrated out of trajectories to save space.

function meshbook()
% Create Mesh from 2D shear data
    
    clear MeshBook MeshList
    clear book4
    
    book4 = book1;
    book4(~book4) = NaN;
    
    MeshBook = book4(1:2,1:25,1:121);
    MeshBook(2,:,:) = MeshBook(2,:,:)*-1;
    node = 0;
    
    for i = 1:size(MeshBook,2)
        for j = 1:size(MeshBook,3)
            node = node + 1;
            if isnan(MeshBook(1,i,j)) == 1
                MeshBook(1,i,j) = MeshBook(1,i-1,j);
            end
            if isnan(MeshBook(2,i,j)) == 1
                MeshBook(2,i,j) = MeshBook(2,i-1,j);
            end
            MeshBook(3,i,j) = i;
            MeshBook(4,i,j) = node;
            
            
            MeshList(2,node) = MeshBook(1,i,j);
            MeshList(3,node) = MeshBook(2,i,j);
            MeshList(4,node) = MeshBook(3,i,j);
            MeshList(1,node) = MeshBook(4,i,j);
            
            MeshList2(2,node) = MeshBook(1,1,j);
            MeshList2(3,node) = MeshBook(2,1,j);
            MeshList2(4,node) = MeshBook(3,i,j);
            MeshList2(1,node) = MeshBook(4,i,j);
        end
    end
    
    %Building Elements
    %Start with the first frame
    clear EleList EleList2
    Element = 0;
    for i = 1:size(MeshBook,3)
        clear current diff
        current(1:2,:) = MeshBook(1:2,1,:);
        diff(1,:) = current(1,:) - MeshBook(1,1,i);
        diff(2,:) = current(2,:) - MeshBook(2,1,i);
        diff(3,:) = sqrt(diff(1,:).^2.+diff(2,:).^2);
        diff(diff<-5) = NaN;
        diff(diff>22) = NaN;
        for j = 1:size(diff,2)
            if isnan(diff(1,j))==1 || isnan(diff(2,j))== 1 || isnan(diff(3,j)) ==1 || diff(3,j) == 0
                diff(1:3,j) = 0;
            end
        end
        
        clear current2
        current2 = find(diff(3,:));
        
        if size(current2,2) == 3
            current3 = diff(1:3,current2);
            mesh3 = current2(find(current3(3,:)==max(current3(3,:))));
            mesh4 = current2(find(current3(2,:)==min(current3(2,:))));
            mesh2 = current2(find(current3(1,:)==min(current3(1,:))));
            for k = 1:(size(MeshBook,2)-1)
                Element = Element + 1;
                
                EleList(i,1,k) = MeshBook(4,k,i); %5
                EleList(i,2,k) = MeshBook(4,k,mesh2); %6
                EleList(i,3,k) = MeshBook(4,k,mesh3); %7
                EleList(i,4,k) = MeshBook(4,k,mesh4); %8
                EleList(i,5,k) = MeshBook(4,k+1,i); %1
                EleList(i,6,k) = MeshBook(4,k+1,mesh2); %2
                EleList(i,7,k) = MeshBook(4,k+1,mesh3); %3
                EleList(i,8,k) = MeshBook(4,k+1,mesh4); %4
                
                EleList(i,9,k) = Element;
                EleList2(2:9,Element) = EleList(i,1:8,k);
                EleList2(1,Element) = Element;
            end
        end
        
        
    end
    
    %Write the deformed mesh
    
    meshTxt = fopen('mesh.txt','wt');
    nodesFormat = '<node id=" %d "> %f, %f, %f </node>\n';
    fprintf(meshTxt,'<?xml version="1.0" encoding="ISO-8859-1"?>\n<febio_spec version="2.5">\n<Geometry>\n<Nodes name="Part1">\n');
    fprintf(meshTxt,nodesFormat,MeshList(1:4,:));
    fprintf(meshTxt,'</Nodes>\n<Elements type="hex8" mat="1" name="part1">\n');
    elementsFormat = '<elem id=" %d "> %d, %d, %d, %d, %d, %d, %d, %d </elem>\n';
    fprintf(meshTxt,elementsFormat,EleList2(1:9,:));
    fprintf(meshTxt,'</Elements>\n</Geometry></febio_spec>');
    fclose(meshTxt);
    
    %Write the undeformed Mesh
    
    meshTxt = fopen('meshUD.txt','wt');
    nodesFormat = '<node id=" %d "> %f, %f, %f </node>\n';
    fprintf(meshTxt,'<?xml version="1.0" encoding="ISO-8859-1"?>\n<febio_spec version="2.5">\n<Geometry>\n<Nodes name="Part1">\n');
    fprintf(meshTxt,nodesFormat,MeshList2(1:4,:));
    fprintf(meshTxt,'</Nodes>\n<Elements type="hex8" mat="1" name="part1">\n');
    elementsFormat = '<elem id=" %d "> %d, %d, %d, %d, %d, %d, %d, %d </elem>\n';
    fprintf(meshTxt,elementsFormat,EleList2(1:9,:));
    fprintf(meshTxt,'</Elements>\n</Geometry></febio_spec>');
    fclose(meshTxt);
    
    %Display the nodes from the deformed mesh
    
    figure
    hold on
    for i = 1:size(MeshBook,2)
        scatter3(MeshBook(1,i,:),MeshBook(2,i,:),i*ones(1,size(MeshBook,3)),'.','g');
        scatter3(MeshBook(1,1,:),MeshBook(2,1,:),i*ones(1,size(MeshBook,3)),'.','b');
    end
    hold off
    
end
end