import gzip, re
import pandas as pd


def list_to_dict(l):
    """Convert list to dict."""
    return {k: v for k, v in (x.split('=') for x in l)}


if __name__ == "__main__":

    # store column info
    cv_columns = {}
    with gzip.open('clinvar.vcf.gz', 'rt') as f:
        for metaline in f:
            if metaline.startswith('##INFO'):
                colname = re.search('ID=(\w+),',
                                    metaline.strip('#\n'))
                coldesc = re.search('.*Description=(.*)>',
                                    metaline.strip('#\n'))
                cv_columns[colname.group(1)] = coldesc.group(1).strip('"')

    # read clinvar vcf
    cv_df = pd.read_csv('clinvar.vcf.gz', sep='\t', comment='#', header=None,
                        usecols=[0, 1, 2, 3, 4, 7], dtype={0: object})

    # convert dictionaries to columns
    cv_df = pd.concat([cv_df.drop([7], axis=1),
                       cv_df[7].str.split(';')
                            .apply(list_to_dict)
                            .apply(pd.Series)], axis=1
                      )
    # rename columns
    cv_df.rename(columns={0: 'CHROM',
                          1: 'POS',
                          2: 'ID',
                          3: 'REF',
                          4: 'ALT'},
                 inplace=True)

    # assign classes
    cv_df['CLASS'] = 0
    cv_df.set_value(cv_df.CLNSIGCONF.notnull(), 'CLASS', 1)

    # convert NaN to 0 where allele frequencies are null
    cv_df[['AF_ESP', 'AF_EXAC', 'AF_TGP']] = \
        cv_df[['AF_ESP', 'AF_EXAC', 'AF_TGP']].fillna(0)

    # select variants that have beeen submitted by multiple organizations.
    cv_df = \
        cv_df.loc[cv_df['CLNREVSTAT']
                  .isin(['criteria_provided,_multiple_submitters,_no_conflicts',
                        'criteria_provided,_conflicting_interpretations'])]

    # Reduce the size of the dataset below
    cv_df.drop(columns=['ALLELEID', 'RS', 'DBVARID'], inplace=True)
    # drop columns that would reveal class
    cv_df.drop(columns=['CLNSIG', 'CLNSIGCONF', 'CLNREVSTAT'], inplace=True)
    # drop these redundant columns
    cv_df.drop(columns=['ID', 'CLNVCSO'], inplace=True)

    cv_df.to_csv('clinvar_conflicting.csv', index=False)
