process SanitizeProteins {
  tag "sanitize"
  cpus 1
  input:
    path proteins
  output:
    path "proteins.clean.faa", emit: cleaned

  """
  awk '
    BEGIN{
      valid="ACDEFGHIKLMNPQRSTVWYBXZOU";
      split(valid,ok,"");
      for(i in ok) map[ok[i]]=1
    }
    /^>/{
      if (hdr != "" && seqlen >= 10) { print hdr; print seq }
      hdr=\$0; seq=""; seqlen=0; next
    }
    {
      gsub(/[^A-Za-z]/,"")
      t=toupper(\$0)
      out=""
      for (i=1; i<=length(t); i++){
        c=substr(t,i,1)
        if (map[c]) out=out c
      }
      seq=seq out
      seqlen=length(seq)
    }
    END{
      if (hdr != "" && seqlen >= 10) { print hdr; print seq }
    }
  ' "${proteins}" > proteins.clean.faa

  test -s proteins.clean.faa
  """
}
