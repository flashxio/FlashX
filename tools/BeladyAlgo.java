import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.PriorityQueue;
import java.util.Comparator;
import java.util.Vector;

public class BeladyAlgo {
	private int shift = 0;
	private int cacheSize;		// in pages

	private HashMap<Integer, Access> accessMap = new HashMap<Integer, Access>();
	private PriorityQueue<Access> accessQueue;
	
	public void printCache() {
		Vector<Access> accesses = new Vector<Access>();
		for (Iterator<Access> it = accessQueue.iterator(); it.hasNext(); )
			accesses.add(it.next());
		Collections.sort(accesses, new Comparator<Access>() {
			@Override
			public int compare(Access a0, Access a1) {
				return (int) (a0.off - a1.off);
			}
		});
		for (int i = 0; i < accesses.size(); i++) {
			System.out.print(accesses.get(i) + " ");
		}
		System.out.println();
	}
	
	public BeladyAlgo(int cacheSize) {
		this.cacheSize = cacheSize;
		accessQueue = new PriorityQueue<Access>(cacheSize, new Comparator<Access>(){
			@Override
			public int compare(Access a1, Access a2) {
				return a2.seq - a1.seq;
			}
			
		});
	}
	
	private class Access {
		public int seq;
		private int off;
		public Access(int off, int seq) {
			this.off = off;
			this.seq = seq;
		}
		
		@Override
		public boolean equals(Object obj) {
			Access access = (Access) obj;
			return access.off == this.off;
		}
		
		@Override
		public String toString() {
			return "(" + off + "," + seq +")"; 
		}
	}
	
	public static int find(int offs[], int shift, int off) {
		for (int i = shift; i < offs.length; i++) {
			if (offs[i] == off)
				return i;
		}
		return offs.length;
	}
	
	public void cacheHit(int off, int offs[]) {
		Access access = accessMap.get(off);
		access.seq = find(offs, shift + 1, off);
		accessQueue.remove(access);
		accessQueue.offer(access);
	}
	
	public void cacheMiss(int off, int offs[]) {
		Access access = accessQueue.peek();
		int nextIdx = find(offs, shift + 1, off);
		if (nextIdx > access.seq)
			return;
		accessQueue.poll();
		accessMap.remove(access.off);
		access = new Access(off, nextIdx);
		accessQueue.offer(access);
		accessMap.put(off, access);
	}
	
	public void addToCache(int off, int offs[]) {
		int nextIdx = find(offs, shift + 1, off);
		Access access = new Access(off, nextIdx);
		accessQueue.offer(access);
		accessMap.put(off, access);
	}
	
	public int access(int offs[]) {
		int cacheHits = 0;
		while (shift < offs.length) {
			int off = offs[shift];
			if (accessMap.containsKey(off)) {
				cacheHit(off, offs);
				cacheHits++;
			}
			else if (accessMap.size() < cacheSize) {
				addToCache(off, offs);
			}
			else {
				cacheMiss(off, offs);
			}
			System.out.println("map size: " + accessMap.size());
			shift++;
		}
		return cacheHits;
	}
	
	public static int[] loadWorkloadOffset(String filename) throws IOException {
		File file = new File(filename);
		long size = file.length();
		int offs[] = new int[((int) size) / 8];
		byte bytes[] = new byte[8];
		byte sizeBytes[] = new byte[4];
		byte padding[] = new byte[4];
		BufferedInputStream in = new BufferedInputStream(new FileInputStream(filename));
		int idx = 0;
		while (in.available() > 0) {
			in.read(bytes);
			in.read(sizeBytes);
			in.read(padding);
			long longValue = (long)(bytes[7] & 0xFF) << 56
		         | (long)(bytes[6] & 0xFF) << 48
		         | (long)(bytes[5] & 0xFF) << 40
		         | (long)(bytes[4] & 0xFF) << 32
		         | (long)(bytes[3] & 0xFF) << 24
		         | (long)(bytes[2] & 0xFF) << 16
		         | (long)(bytes[1] & 0xFF) <<  8
		         | (long)(bytes[0] & 0xFF) <<  0;
			int length = (int)(sizeBytes[3] & 0x7F) << 24
		         | (int)(sizeBytes[2] & 0xFF) << 16
		         | (int)(sizeBytes[1] & 0xFF) <<  8
		         | (int)(sizeBytes[0] & 0xFF) <<  0;
			long end = longValue + length;
			while (longValue <= end) {
				offs[idx] = (int) (longValue / 4096);
				longValue = ((long) offs[idx]) * 4096 + 4096;
				idx++;
			}
		}
		in.close();
		int tmp[] = new int[idx];
		System.arraycopy(offs, 0, tmp, 0, idx);
		return tmp;
	}
	
	public static void main(String args[]) throws IOException {
		if (args.length < 2) {
			System.err.println("BeladyAlgo file cacheSize");
			return;
		}
		int offs[] = loadWorkloadOffset(args[0]);
		System.out.println("load all long integers");
		int cacheSize = Integer.parseInt(args[1]);		// in number of pages.
		BeladyAlgo algo = new BeladyAlgo(cacheSize);
		int hits = algo.access(offs);
		System.out.println("There are " + hits + " out of " + offs.length + " accesses");
	}
}
