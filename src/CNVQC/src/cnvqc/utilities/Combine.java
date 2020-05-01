package cnvqc.utilities;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

//Small class that functions like perl's 'join()' method
public class Combine {

	public static String combineList(List<String> s, String glue) {
	  int k = s.size();
	  if ( k == 0 )
	  {
	    return null;
	  }
	  StringBuilder out = new StringBuilder();
	  out.append( s.get(0) );
	  for ( int x=1; x < k; ++x )
	  {
	    out.append(glue).append(s.get(x));
	  }
	  return out.toString();
	}
	public static String combineArray(String[] s, String glue) {
	  int k = s.length;
	  if ( k == 0 ) {
	    return null;
	  }
	  StringBuilder out = new StringBuilder();
	  out.append( s[0] );
	  for ( int x=1; x < k; ++x ) {
	    out.append(glue).append(s[x]);
	  }
	  return out.toString();
	}
	public static String combineSet(Set<String> s, String glue) {
		int k = s.size();
		if ( k == 0 )
		{
			return null;
		}
		List<String> t = new ArrayList<String>();
		for (String str : s) {
			t.add(str);
		}
		return combineList(t, glue);
	}
}
