package PrescribedModularity;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import delft.recursiveExploration.Basics.GraphSparse;

public class modularityApproachingMixtureOneTwo {
	
	public int L; // given number of links
	public int c; // given number of partitions
	public double M; // given modularity
	public double toler; // tolerance
	
	int min_i;
    int min_j;
    
    int L_int;
	
	public int[][] C; // community graph
	public int[] internalLinks; //internal links
	public GraphSparse G; // the graph
	
	public int MinDistance(int [] ar){
        int [] a = ar;
        int aSize = a.length;
        int dMin = Integer.MAX_VALUE;//MaxInt
        min_i = 1;
        min_j = 2;
        
        for(int i=1; i<aSize; i++)
        {
            for(int j=i+1; j<aSize;j++)
            {
            	if((a[i]-a[j])>2)
            	{
            		if((a[i]-a[j])<dMin)
            		{
            			dMin = a[i]-a[j];
            			min_i = i;
            		    min_j = j;
            		}
            	}
            	else if((a[j]-a[i])>2)
            	{
            		if((a[j]-a[i])<dMin)
            		{
            			dMin = a[j]-a[i];
            			min_i = j;
            		    min_j = i;
            		}
            	}           	
            }
        }
        
        return dMin;
    }
	
	public modularityApproachingMixtureOneTwo(int _L, int _c, double _M, double _t) throws IOException
	{
		L = _L;
		c = _c;
		M = _M;
		toler = _t;
		
		System.out.println("------------------------ Random --------------------------");
		System.out.println();
		
		double max_m = initialize();		
		if(approxAlgorithm(max_m))
		{
			System.out.println("There is a graph with " + L + " links, partitioned into " + c + " communities and modularity " + M);
			System.out.println();
			
			for(int i = 1; i<=c; i++)
				System.out.print(i +"comm\t");
			
			System.out.println("L_int (inter-communities links)");
			
			for(int i = 1; i<=c; i++)
				System.out.print(internalLinks[i] +"\t");
			
			System.out.println(L_int);
			
			System.out.println();
			System.out.println("Inter-community links:");
			
			for(int k = 1; k<=c; k++)
				for(int r = k+1; r<=c; r++)
					if(C[k][r]>0)
					{
						System.out.println(k+" comm. -> " + r +" comm. " + C[k][r] + " link(s)" +"\t");
					}
					
			
		}
		else
		{
			System.out.println("There is no graph with " + L + " links, partitioned into " + c + " communities and modularity " + M);
		}
		
		System.out.println();
	}


	private double initialize() {
		
		C = new int[c+1][c+1];
		internalLinks = new int[c+1];

		for(int k=1; k<=c; k++)
		{
			internalLinks[k] = 0;
			
			for(int r=k; r<=c; r++)
			{
				C[r][k] = 0;
				C[k][r] = 0;
			}
		}
		
		int r = L%c;
		int k = L/c;
		double mmax = 1.0 - 1.0/c + ((double)(c-1))/L;
		
		internalLinks[1] = k; int i = 2;
		
		if(r==0)
		{
			while(i<=c)
			{
				C[i-1][i] = C[i][i-1] = 1;
				internalLinks[i] = k-1;
				i++;
			}
			
			mmax -= 1.0/(2*L*L);
		}
		else if(r<(c/2))
		{
			while(i<=c-r)
			{
				C[i-1][i] = C[i][i-1] = 1;
				internalLinks[i] = k-1;
				
				if(i<=r)
				{
					C[c-i+1][i] = C[i][c-i+1] = 1;
					internalLinks[c-i+1] = k;
				}
				
				i++;
			}
			
			internalLinks[i] = k;
			mmax -= ((double)(r*(c-2*r)))/(2*c*L*L);
		}
		else
		{
			while(i<=r)
			{
				C[i-1][i] = C[i][i-1] = 1;
				internalLinks[i] = k;
				
				if(i<=c-r)
				{
					C[c-i+1][i] = C[i][c-i+1] = 1;
					internalLinks[i]--;
					internalLinks[c-i+1] = k;
				}
				
				i++;
			}
			
			mmax -= ((double)((c-r)*(2*r-c)))/(2*c*L*L);
		}
		
		return mmax;
	}
	
	private boolean approxAlgorithm(double max_m) throws IOException {
		
		double curr_m = max_m;
		
		if(curr_m-toler<M)
			return false;
		

		
		int cur_i = -1;
		int cur_j = -1;
		int state = 0;
		
		int[] cumDeg = new int[c+1];
		
		for(int k=1; k<=c; k++)
		{
			cumDeg[k] = 2*internalLinks[k];
			
			for(int r=1; r<=c; r++)
				cumDeg[k]+=C[k][r];
		}
			
		
		//phase 2		
		while(Math.abs(curr_m-M)>toler)
		{			
			MinDistance(cumDeg);
			
			double rand = Math.random();
			
			if(rand<0.5)
			{
				if(curr_m>M)
				{
					if((state==2)&&(min_i==cur_i)&&(min_j==cur_j))
						return false;
					
					if(internalLinks[min_j]<20)
						break;
					
					curr_m-=delta_2(min_i,min_j);
					
					internalLinks[min_i]++; cumDeg[min_i]+=2;
					internalLinks[min_j]--; cumDeg[min_j]-=2;
					state = 1;
					
					//System.out.println("m = " + curr_m + " decr; phase 1");
				}
				else
				{
					if((state==1)&&(min_i==cur_i)&&(min_j==cur_j))
						return false;
					
					if(internalLinks[min_i]<20)
						break;
					
					curr_m+=delta_2(min_i,min_j);
					
					internalLinks[min_i]--; cumDeg[min_i]-=2;
					internalLinks[min_j]++; cumDeg[min_j]+=2;
					state = 2;
					
					//System.out.println("m = "+curr_m + " incr; phase 1");
				}
				
				cur_i = min_i;
				cur_j = min_j;
			}
			else
			{
				
				if(curr_m>M)
				{
					if((state==2)&&(min_i==cur_i)&&(min_j==cur_j))
						return false;
					
					if(internalLinks[min_j]<20)
						break;
					
					curr_m-=delta_1(min_i,min_j);
					C[min_i][min_j]++;
					C[min_j][min_i]++;
					internalLinks[min_j]--;
					
					cumDeg[min_i]++;
					cumDeg[min_j]--;
					
					state = 1;
					
					//System.out.println("m = " + curr_m + " decr; phase 2");
				}
				else
				{
					if((state==1)&&(min_i==cur_i)&&(min_j==cur_j))
						return false;
					
					curr_m+=delta_1(min_i,min_j);
					C[min_i][min_j]--;
					C[min_j][min_i]--;
					internalLinks[min_j]++;
					
					cumDeg[min_i]--;
					cumDeg[min_j]++;
					
					state = 2;
					
					//System.out.println("m = " + curr_m + " incr; phase 2");
				}
				
				cur_i = min_i;
				cur_j = min_j;
				
			}		
		}
				
		L_int = 0;
		
		for(int k=1; k<=c; k++)
		{	
			for(int r=k+1; r<=c; r++)
					L_int += C[k][r];
		}
				
		return true;		
	}

	private double delta_1(int i, int j) {
		
		double d_i = 2*internalLinks[i];
		double d_j = 2*internalLinks[j];
		
		for(int k=1; k<=c; k++)
		{
			d_i += (double)(C[k][i]);
			d_j += (double)(C[k][j]);
		}
		
		return ((double)(2*L+d_j-d_i-1))/(2*L*L);
	}

	private double delta_2(int i, int j) {
		
		double d_i = 2*internalLinks[i];
		double d_j = 2*internalLinks[j];
		
		for(int k=1; k<=c; k++)
		{
			d_i += (double)(C[k][i]);
			d_j += (double)(C[k][j]);
		}
		
		return ((double)(d_i-d_j-2))/(L*L);
	}

}
