#' Create a table of types of point mutation
#' A type of point mutation is TAA> TGA (triplet = "TGA", mutation ="A" )
#' @param triplets E.g. TCG, ACT, AAA etc. These might be recovered using the expandSeq function on a VCF file.
#' @param mutations The mutation at the centre of the triplet e.g. TCG > TGG is mutation G
#' @keywords Poisson, mutants, kataegesis
#' @export
#' @examples
#' tripletmutcount(VCF, triplets = triplets, mutations = MUT)
#'

tripletmut.count= function(VCF, triplets = triplets, mutations = MUT)

{
  triplets = deparse(substitute(triplets))
  mutations = deparse(substitute(mutations))
  
  #prep
  all.mutations <- c("A","C", "G", "T")
  all.triplets <- paste0(rep(all.mutations, times=16), rep(all.mutations, each=4), rep(all.mutations, each=16))

  # This is awkward but I need to make a table of all possible mutation types and merge with the actual mutation table
  # remember a mutaion canot be AAA>A or CTT>T - so these are removed now.
  mutation.allposs <- as.tbl(data.frame(triplets=rep(all.triplets,each=4), mutations=all.mutations)) %>%
    filter(substr(triplets,2,2) != mutations) %>%
    arrange(triplets, mutations)

  # this is the actual mutation type table
  mutation.table <- table(paste0(VCF$triplets, VCF$mutations))
  mutation.actual = as.tbl(data.frame(triplets=substr(names(mutation.table),1,3), mutations=substr(names(mutation.table),4,4), count=as.vector(mutation.table)))

  #VCF %>% group_by(triplets, mutations)  %>% summarise_each(funs(length), mutations)

  # this is just in case there are missing mutations in the actual data, to get 0 counts for them rather than blank rows.
  mutation.join = left_join(mutation.allposs, mutation.actual) %>%
    mutate(count = ifelse(is.na(count), 0, count))

  # each mutation has a single reverse complement which is identical from a molecular biol. viewpoint
  # so really there are 96 unique point mutation NOT 192
  mutation.join <- mutation.join %>%
    mutate(rctriplets = as.character(reverseComplement(DNAStringSet(triplets)))) %>%
    mutate(rcmutations = as.character(reverseComplement(DNAStringSet(mutations))))

  # this code adds the counts from the reverse complement then removes the 96 redundant rows.
  matches <- with(mutation.join, match(paste0(triplets,mutations), paste0(rctriplets,rcmutations)))
  mutation.join <- mutate(mutation.join, count.both= count + count[matches])
  matches <- matches[matches < 1: length(matches)]

  # tidy the data
  mutation.join <- as.tbl(data.frame(mutation.join[matches,])) %>%
    arrange(triplets,mutations)

  return(mutation.join)
}
