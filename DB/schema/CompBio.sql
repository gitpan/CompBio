# MySQL dump 8.16
#
# Host: sigler    Database: CompBio
#--------------------------------------------------------
# Server version	3.23.41

#
# Table structure for table 'Correlation'
#

DROP TABLE IF EXISTS Correlation;
CREATE TABLE Correlation (
  Genome INT UNSIGNED zerofill default NULL,
  Swissprot INT unsigned zerofill default NULL,
  Genbank INT unsigned zerofill default NULL,
  PDB INT unsigned zerofill default NULL,
  KEY Genome (Genome),
  KEY Swissprot (Swissprot),
  KEY Genbank (Genbank),
  KEY PDB (PDB)
) TYPE=MyISAM PACK_KEYS=1 CHECKSUM=1 COMMENT='All DBKeys know for a specific gene object';

#
# Table structure for table 'ID'
#
# Suffix is attached to the name value, such as version # for genbanck accesion #
# or chain id for PDB

DROP TABLE IF EXISTS ID;
CREATE TABLE ID (
  DBKey INT unsigned zerofill NOT NULL AUTO_INCREMENT,
  Name varchar(25) binary NOT NULL default '',
  Suffix VARCHAR(2),
  DB enum('Genome','Swissprot','Genbank','PDB','EST') NOT NULL,
  TS timestamp(10) NOT NULL,
  AlternateNames varchar(50) default NULL,
  PRIMARY KEY  (DBKey),
  KEY Name (Name)
) TYPE=MyISAM PACK_KEYS=1 CHECKSUM=1 COMMENT='Generate db key for gene objects. AlternateNames are csv''s';

###
### MERGE TABLES - edit paths by hand after creation
### Need to create temp table for merge

# blast tables

DROP TABLE IF EXISTS gBlast;
CREATE TABLE gBlast (
  ID INT UNSIGNED ZEROFILL NOT NULL,
  Match_ID INT UNSIGNED ZEROFILL NOT NULL,
  EffectiveIdentities float(8,3) unsigned NOT NULL default '0.000',
  EquivalentIdentities tinyint(3) unsigned NOT NULL default '0',
  Qstart smallint(5) unsigned default NULL,
  Qend smallint(5) unsigned default NULL,
  Mstart smallint(5) unsigned default NULL,
  Mend smallint(5) unsigned default NULL,
  PRIMARY KEY (ID,Match_ID)
) TYPE=MyISAM PACK_KEYS=1 CHECKSUM=1 COMMENT="Corrilation using Blast";

DROP TABLE IF EXISTS sBlast;
CREATE TABLE sBlast (
  ID INT UNSIGNED ZEROFILL NOT NULL,
  Match_ID INT UNSIGNED ZEROFILL NOT NULL,
  EffectiveIdentities float(8,3) unsigned NOT NULL default '0.000',
  EquivalentIdentities tinyint(3) unsigned NOT NULL default '0',
  Qstart smallint(5) unsigned default NULL,
  Qend smallint(5) unsigned default NULL,
  Mstart smallint(5) unsigned default NULL,
  Mend smallint(5) unsigned default NULL,
  PRIMARY KEY (ID,Match_ID)
) TYPE=MyISAM PACK_KEYS=1 CHECKSUM=1 COMMENT="Corrilation using Blast";

DROP TABLE IF EXISTS nBlast;
CREATE TABLE nBlast (
  ID INT UNSIGNED ZEROFILL NOT NULL,
  Match_ID INT UNSIGNED ZEROFILL NOT NULL,
  EffectiveIdentities float(8,3) unsigned NOT NULL default '0.000',
  EquivalentIdentities tinyint(3) unsigned NOT NULL default '0',
  Qstart smallint(5) unsigned default NULL,
  Qend smallint(5) unsigned default NULL,
  Mstart smallint(5) unsigned default NULL,
  Mend smallint(5) unsigned default NULL,
  PRIMARY KEY (ID,Match_ID)
) TYPE=MyISAM PACK_KEYS=1 CHECKSUM=1 COMMENT="Corrilation using Blast";

DROP TABLE IF EXISTS pBlast;
CREATE TABLE pBlast (
  ID INT UNSIGNED ZEROFILL NOT NULL,
  Match_ID INT UNSIGNED ZEROFILL NOT NULL,
  EffectiveIdentities float(8,3) unsigned NOT NULL default '0.000',
  EquivalentIdentities tinyint(3) unsigned NOT NULL default '0',
  Qstart smallint(5) unsigned default NULL,
  Qend smallint(5) unsigned default NULL,
  Mstart smallint(5) unsigned default NULL,
  Mend smallint(5) unsigned default NULL,
  PRIMARY KEY (ID,Match_ID)
) TYPE=MyISAM PACK_KEYS=1 CHECKSUM=1 COMMENT="Corrilation using Blast";

DROP TABLE IF EXISTS Blast;
CREATE TABLE Blast (
  ID INT UNSIGNED ZEROFILL NOT NULL,
  Match_ID INT UNSIGNED ZEROFILL NOT NULL,
  EffectiveIdentities FLOAT(8,3) UNSIGNED NOT NULL,
  EquivalentIdentities TINYINT UNSIGNED NOT NULL,
  Qstart SMALLINT UNSIGNED DEFAULT NULL,
  Qend SMALLINT UNSIGNED DEFAULT NULL,
  Mstart SMALLINT UNSIGNED DEFAULT NULL,
  Mend SMALLINT UNSIGNED DEFAULT NULL,
  KEY (ID,Match_ID)
) TYPE=MERGE UNION=(sBlast,gBlast,nBlast,pBlast) COMMENT="Corrilation using Blast";

DROP TABLE IF EXISTS sBlast;
DROP TABLE IF EXISTS gBlast;
DROP TABLE IF EXISTS nBlast;
DROP TABLE IF EXISTS pBlast;

# annotation tables

DROP TABLE IF EXISTS sAnnotation;
CREATE TABLE sAnnotation (
  ID INT UNSIGNED ZEROFILL NOT NULL,
  Definition text NOT NULL,
  Supplemental text,
  Submitted text,
  PRIMARY KEY  (ID)
) TYPE=MyISAM PACK_KEYS=1 CHECKSUM=1 COMMENT="Original and extra annotation";

DROP TABLE IF EXISTS gAnnotation;
CREATE TABLE gAnnotation (
  ID INT UNSIGNED ZEROFILL NOT NULL,
  Definition text NOT NULL,
  Supplemental text,
  Submitted text,
  PRIMARY KEY  (ID)
) TYPE=MyISAM PACK_KEYS=1 CHECKSUM=1 COMMENT="Original and extra annotation";

DROP TABLE IF EXISTS nAnnotation;
CREATE TABLE nAnnotation (
  ID INT UNSIGNED ZEROFILL NOT NULL,
  Definition text NOT NULL,
  Supplemental text,
  Submitted text,
  PRIMARY KEY  (ID)
) TYPE=MyISAM PACK_KEYS=1 CHECKSUM=1 COMMENT="Original and extra annotation";

DROP TABLE IF EXISTS pAnnotation;
CREATE TABLE pAnnotation (
  ID INT UNSIGNED ZEROFILL NOT NULL,
  Definition text NOT NULL,
  Supplemental text,
  Submitted text,
  PRIMARY KEY  (ID)
) TYPE=MyISAM PACK_KEYS=1 CHECKSUM=1 COMMENT="Original and extra annotation";

DROP TABLE IF EXISTS Annotation;
CREATE TABLE Annotation (
  ID INT UNSIGNED ZEROFILL NOT NULL,
  Definition TEXT NOT NULL,
  Supplemental TEXT,
  Submitted TEXT,
  PRIMARY KEY  (ID)
) TYPE=MERGE UNION=(sAnnotation,gAnnotation,nAnnotation,pAnnotation) COMMENT="Original and extra annotation";

DROP TABLE IF EXISTS sAnnotation;
DROP TABLE IF EXISTS gAnnotation;
DROP TABLE IF EXISTS nAnnotation;
DROP TABLE IF EXISTS pAnnotation;

# peptide tables

DROP TABLE IF EXISTS sPeptide;
CREATE TABLE sPeptide (
  ID INT UNSIGNED ZEROFILL NOT NULL,
  Sequence text NOT NULL,
  PRIMARY KEY  (ID)
) TYPE=MyISAM PACK_KEYS=1 CHECKSUM=1;

DROP TABLE IF EXISTS gPeptide;
CREATE TABLE gPeptide (
  ID INT UNSIGNED ZEROFILL NOT NULL,
  Sequence text NOT NULL,
  PRIMARY KEY  (ID)
) TYPE=MyISAM PACK_KEYS=1 CHECKSUM=1;

DROP TABLE IF EXISTS nPeptide;
CREATE TABLE nPeptide (
  ID INT UNSIGNED ZEROFILL NOT NULL,
  Sequence text NOT NULL,
  PRIMARY KEY  (ID)
) TYPE=MyISAM PACK_KEYS=1 CHECKSUM=1;

DROP TABLE IF EXISTS pPeptide;
CREATE TABLE pPeptide (
  ID INT UNSIGNED ZEROFILL NOT NULL,
  Sequence text NOT NULL,
  PRIMARY KEY  (ID)
) TYPE=MyISAM PACK_KEYS=1 CHECKSUM=1;

DROP TABLE IF EXISTS Peptide;
CREATE TABLE Peptide (
  ID INT UNSIGNED ZEROFILL NOT NULL,
  Sequence TEXT NOT NULL,
  PRIMARY KEY  (ID)
) TYPE=MERGE UNION=(sPeptide,gPeptide,nPeptide,pPeptide) PACK_KEYS=1 CHECKSUM=1;

DROP TABLE IF EXISTS sPeptide;
DROP TABLE IF EXISTS gPeptide;
DROP TABLE IF EXISTS nPeptide;
DROP TABLE IF EXISTS pPeptide;

# keyword tables

DROP TABLE IF EXISTS sKeywords;
CREATE TABLE sKeywords (
  ID INT UNSIGNED ZEROFILL NOT NULL,
  GeneName varchar(12) default NULL,
  Function varchar(32) default NULL,
  PRIMARY KEY  (ID),
  KEY (GeneName)
) TYPE=MyISAM PACK_KEYS=1 CHECKSUM=1 COMMENT="Single term fields, such as HisH or amidotransferase";

DROP TABLE IF EXISTS gKeywords;
CREATE TABLE gKeywords (
  ID INT UNSIGNED ZEROFILL NOT NULL,
  GeneName varchar(12) default NULL,
  Function varchar(32) default NULL,
  PRIMARY KEY  (ID),
  KEY (GeneName)
) TYPE=MyISAM PACK_KEYS=1 CHECKSUM=1 COMMENT="Single term fields, such as HisH or amidotransferase";

DROP TABLE IF EXISTS nKeywords;
CREATE TABLE nKeywords (
  ID INT UNSIGNED ZEROFILL NOT NULL,
  GeneName varchar(12) default NULL,
  Function varchar(32) default NULL,
  PRIMARY KEY  (ID),
  KEY (GeneName)
) TYPE=MyISAM PACK_KEYS=1 CHECKSUM=1 COMMENT="Single term fields, such as HisH or amidotransferase";

DROP TABLE IF EXISTS pKeywords;
CREATE TABLE pKeywords (
  ID INT UNSIGNED ZEROFILL NOT NULL,
  GeneName varchar(12) default NULL,
  Function varchar(32) default NULL,
  PRIMARY KEY  (ID),
  KEY (GeneName)
) TYPE=MyISAM PACK_KEYS=1 CHECKSUM=1 COMMENT="Single term fields, such as HisH or amidotransferase";

DROP TABLE IF EXISTS Keywords;
CREATE TABLE Keywords (
  ID INT UNSIGNED ZEROFILL NOT NULL,
  GeneName VARCHAR(12) DEFAULT NULL,
  Function VARCHAR(32) DEFAULT NULL,
  PRIMARY KEY  (ID),
  KEY (GeneName)
) TYPE=MERGE UNION=(sKeywords,gKeywords,nKeywords,pKeywords) PACK_KEYS=1 CHECKSUM=1 COMMENT="Single term fields, such as HosH or amidotransferase";

DROP TABLE IF EXISTS sKeywords;
DROP TABLE IF EXISTS gKeywords;
DROP TABLE IF EXISTS nKeywords;
DROP TABLE IF EXISTS pKeywords;

