function resetVertex ( ax )
  % extract the x axis vertext data
  % X, Y and Z row of the start and end of the individual axle.
  ax.XAxis.Axle.VertexData(1,1) = 0;
  % repeat for Y (set 2nd row)
  ax.YAxis.Axle.VertexData(2,1) = 0;

  % You can modify the minor Tick values by modifying the vertex data
  % for them, e.g. remove any minor ticks below 0 
  ax.XAxis.MinorTickChild.VertexData(:,ax.XAxis.MinorTickChild.VertexData(1,:)<0) = [];
  ax.YAxis.MinorTickChild.VertexData(:,ax.YAxis.MinorTickChild.VertexData(2,:)<0) = [];
end