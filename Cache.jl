using SHA
using Serialization

export @cache

macro cache(expr)
	hash = bytes2hex(sha256(string(expr)))
	filename = ".cache/$hash"

	return :(
		try
			deserialize($filename)
		catch
			try
				mkdir(".cache")
			catch
			end
			temp = $expr;
			serialize($filename, temp);
			temp
		end
	)
end
