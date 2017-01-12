/*
 * source code for  
 */

package PrescribedModularity;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import delft.recursiveExploration.Basics.GraphSparse;
import delft.recursiveExploration.Basics.Node;
import delft.recursiveExploration.Basics.Link;
import delft.recursiveExploration.Graph.GraphInfeasibleException;
import delft.recursiveExploration.Graph.RandomGraph;
import delft.recursiveExploration.Metrics.Modularity;

public class generateModularNetworks {
	
	public int L; // given number of links
	public int c; // given number of partitions
	public double M; // given modularity
	public double toler; // tolerance
	
	public modularityApproachingOneThenTwo OneThenTwo;
	public modularityApproachingTwoThenOne TwoThenOne;
	public modularityApproachingMixtureOneTwo MixtureOneTwo;
	public GraphSparse G;
	public HashMap<Integer,ArrayList<Node>> nodeInPartition;
	public HashMap<Node,Node> mapNodeId;
	
	public generateModularNetworks(int _L, int _c, double _M, double _t, int algType) throws IOException, GraphInfeasibleException
	{
		L = _L;
		c = _c;
		M = _M;
		toler = _t;
		String fileSuffix;
		
		RandomGraph randomTemplate;
		
		nodeInPartition = new HashMap<Integer, ArrayList<Node>>(c);
		mapNodeId = new HashMap<Node, Node>();
		
		G = new GraphSparse();
		
		switch(algType)
		{
			case 0:
				MixtureOneTwo = new modularityApproachingMixtureOneTwo(L, c, M, toler);
				fileSuffix = "Random";
				break;
			case 1:
				OneThenTwo = new modularityApproachingOneThenTwo(L, c, M, toler);
				fileSuffix = "OneThenTwo";
				break;
			case 2:
				TwoThenOne = new modularityApproachingTwoThenOne(L, c, M, toler);
				fileSuffix = "TwoThenOne";
				break;
			default:
				TwoThenOne = new modularityApproachingTwoThenOne(L, c, M, toler);
				fileSuffix = "any";
				break;
		}
		
		int nodeCounter = 1;
		ArrayList<Node> lstInPart;
		int numLinks = -1;
		int numNodes = -1;
		
		
		for(int par = 1; par<=c; par++)
		{
			switch(algType)
			{
				case 0:
					numLinks = MixtureOneTwo.internalLinks[par];
					numNodes = (int)((numLinks-2*Math.sqrt(2*numLinks)+1)/2);
					randomTemplate = new RandomGraph(numNodes);
					break;
				case 1:
					numLinks = OneThenTwo.internalLinks[par];
					numNodes = (int)((numLinks-2*Math.sqrt(2*numLinks)+1)/2);
					randomTemplate = new RandomGraph(numNodes);
					break;
				case 2:
					numLinks = TwoThenOne.internalLinks[par];
					numNodes = (int)((numLinks-2*Math.sqrt(2*numLinks)+1)/2);
					randomTemplate = new RandomGraph(numNodes);
					break;
				default:
					randomTemplate = new RandomGraph(numNodes);
					break;
			}
			
			randomTemplate.generate_L(numLinks);
			
			lstInPart = new ArrayList<Node>();
			
			for(Node nd : randomTemplate.getAllNodes()) {
				
				Node gNode = G.addNode(String.valueOf(nodeCounter));
				nodeCounter++;
				mapNodeId.put(nd, gNode);
				lstInPart.add(gNode);			
			}
			
			for(Node sd : randomTemplate.getAllNodes()) {
				for(Node dd : randomTemplate.getNeighbors(sd)) {				
					if(randomTemplate.getmapNodetoID().get(sd)<randomTemplate.getmapNodetoID().get(dd))
						G.addLink(mapNodeId.get(sd), mapNodeId.get(dd));
				}
			}
			
			nodeInPartition.put(par,lstInPart);		
		}
		
		switch(algType)
		{
			case 0:
				for(int k = 1; k<MixtureOneTwo.C.length; k++)
					for(int r = (k+1); r<MixtureOneTwo.C.length; r++)
					{
						int m = 0;
						while(m<MixtureOneTwo.C[k][r])
						{
							int a = (int)(Math.random()*nodeInPartition.get(k).size());
							int b = (int)(Math.random()*nodeInPartition.get(r).size());
							
							G.addLink(nodeInPartition.get(k).get(r),nodeInPartition.get(r).get(b));
							m++;
						}
					}
				break;
			case 1:
				for(int k = 1; k<OneThenTwo.C.length; k++)
					for(int r = (k+1); r<OneThenTwo.C.length; r++)
					{
						int m = 0;
						while(m<OneThenTwo.C[k][r])
						{
							int a = (int)(Math.random()*nodeInPartition.get(k).size());
							int b = (int)(Math.random()*nodeInPartition.get(r).size());
							
							G.addLink(nodeInPartition.get(k).get(r),nodeInPartition.get(r).get(b));
							m++;
						}
					}
				break;
			case 2:
				for(int k = 1; k<TwoThenOne.C.length; k++)
					for(int r = (k+1); r<TwoThenOne.C.length; r++)
					{
						int m = 0;
						while(m<TwoThenOne.C[k][r])
						{
							int a = (int)(Math.random()*nodeInPartition.get(k).size());
							int b = (int)(Math.random()*nodeInPartition.get(r).size());
							
							G.addLink(nodeInPartition.get(k).get(a),nodeInPartition.get(r).get(b));
							m++;
						}
					}
				break;
			default:
				for(int k = 1; k<TwoThenOne.C.length; k++)
					for(int r = (k+1); r<TwoThenOne.C.length; r++)
					{
						int m = 0;
						while(m<TwoThenOne.C[k][r])
						{
							G.addLink(nodeInPartition.get(k).get(0),nodeInPartition.get(r).get(0));
							m++;
						}
					}
				break;
		}
		
		G.setName("ClusteredGraph_L"+L+"c"+c+fileSuffix);
		
		saveAdjacencyFile();
		
		DecimalFormat df = new DecimalFormat("#.##");
			
		System.out.println("Modularity quality coefficient (difference from Newman modularity algorithm) K = " + df.format(M/calculateNewmanModularity()*100.00) + "%"); // Modularity quality coefficient - comparison with Newman modularity
		System.out.println("------------------------------------------------------------------");
	}
	
	private void saveAdjacencyFile() throws IOException {
		
		FileWriter fstream;
		BufferedWriter out = null;		

		fstream = new FileWriter("prescribedModularity//"+G.getName()+".txt");
		out = new BufferedWriter(fstream);
		
		for (Link lk : G.getAllLinks())
			out.write(lk.i + "\t" + lk.j + "\n");

		out.close();		
	}
	
	/*
	 * Newman modularity calculation
	 * 
	*/	
	public double calculateNewmanModularity()
	{
		Modularity modularity;
		modularity = new Modularity(G);
		modularity.initialize();
		
		return modularity.getValue();
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 * @throws GraphInfeasibleException 
	 */
	public static void main(String[] args) throws IOException, GraphInfeasibleException {
		
		generateModularNetworks findGraph = new generateModularNetworks(1000, 5, 0.65, 0.01, 2); // arguments: # links; # communities; desired modularity; tolerance; variant {0,1,2}
	}
}
