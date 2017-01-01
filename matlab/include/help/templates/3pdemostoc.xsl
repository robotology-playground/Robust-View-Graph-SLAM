<?xml version="1.0"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
                xmlns:file="java.io.File"
                exclude-result-prefixes="file"
                version="1.0">

  <!-- Copyright 2014-2015 The MathWorks, Inc. -->

  <xsl:output method="html" indent="yes"/>

  <xsl:param name="demosroot"/>
  <xsl:param name="languageDir"/>
  <xsl:param name="helpdir"/>
  <xsl:param name="docpage"/>

  <!-- Create the full product name -->
  <xsl:variable name="productName">
    <xsl:value-of select="/demos/name"/>
    <xsl:choose>
      <xsl:when test="/demos/type='toolbox'"> Toolbox</xsl:when>
      <xsl:when test="/demos/type='blockset'"> Blockset</xsl:when>
    </xsl:choose>
  </xsl:variable>

  <xsl:template match="/demos">
    <ul>
      <li id="demos.xml" class="toc_expanded">
        <div class="node_element">
          <div class="wrapper img_wrapper"> </div>
          <div class="wrapper text_wrapper">
            <a href="{concat($demosroot,'/demos.xml')}" class="corrected_url">
              <xsl:value-of select="$productName"/><xsl:text> Examples (Supplemental Software)</xsl:text>
            </a>
          </div>
        </div>
        <ul>
          <xsl:apply-templates select="child::demosection"/>
        </ul>
      </li>
    </ul>
    <xsl:if test="$helpdir != ''">
      <ul>
        <li id="{$docpage}">
          <div class="node_element">
            <div class="wrapper img_wrapper"> </div>
            <div class="wrapper text_wrapper">
              <a href="{concat($helpdir,'/',$docpage)}" class="corrected_url">Supplemental Software</a>
            </div>
          </div>
        </li>
      </ul>
    </xsl:if>
  </xsl:template>

  <xsl:template match="demosection">
    <!-- g1150253: Counting preceding demo sections to properly render identical sections -->
    <li id="{concat('demos.xml#',translate(label,' ','-'),count(preceding::demosection))}">
      <xsl:if test="demoitem | demosection">
        <xsl:attribute name="class">toc_collapsed</xsl:attribute>
      </xsl:if>
      <div class="node_element">
        <div class="wrapper img_wrapper"> </div>
        <div class="wrapper text_wrapper">
          <a href="{concat($demosroot,'/demos.xml#',translate(label,' ','-'),count(preceding::demosection))}" class="corrected_url">
            <xsl:value-of select="label"/>
          </a>
        </div>
      </div>
      <xsl:if test="demoitem | demosection">
        <ul>
          <xsl:apply-templates select="demoitem | demosection"/>
        </ul>
      </xsl:if>
    </li>
  </xsl:template>

  <xsl:template match="demoitem">
    <!-- Construct title with link -->
    <xsl:variable name="link">
      <xsl:choose>
        <xsl:when test="type and file">
          <xsl:variable name="localFile">
            <xsl:if test="$languageDir != ''">
              <xsl:call-template name="buildLocalizedPath">
                <xsl:with-param name="filepath" select="file"/>
              </xsl:call-template>
            </xsl:if>
          </xsl:variable>
          <xsl:call-template name="checkLocalizedFile">
            <xsl:with-param name="localizedFile" select="$localFile"/>
            <xsl:with-param name="defaultFile" select="file"/>
          </xsl:call-template>
        </xsl:when>
        <xsl:when test="'M-file'=type">
          <xsl:variable name="source" select="source"/>
          <xsl:variable name="localFile">
            <xsl:if test="$languageDir != ''">
              <xsl:value-of select="concat('html/',$languageDir,'/',$source,'.html')"/>
            </xsl:if>
          </xsl:variable>
          <xsl:variable name="defaultFile" select="concat('html/',$source,'.html')"/>
          <xsl:call-template name="checkLocalizedFile">
            <xsl:with-param name="localizedFile" select="$localFile"/>
            <xsl:with-param name="defaultFile" select="$defaultFile"/>
          </xsl:call-template>
        </xsl:when>
        <xsl:otherwise>
          <xsl:text>matlab:indexhelper('</xsl:text>
          <xsl:call-template name="quotify">
            <!-- Strip off the leading file:/// for use with indexhelper. -->
<xsl:with-param name="string" select="substring-after($demosroot,'file:///')"/>
          </xsl:call-template>
          <xsl:text>','</xsl:text>
          <xsl:call-template name="quotify">
            <xsl:with-param name="string" select="source"/>
          </xsl:call-template>
          <xsl:text>','</xsl:text>
          <xsl:call-template name="quotify">
            <xsl:with-param name="string" select="callback"/>
          </xsl:call-template>
          <xsl:text>','</xsl:text>
          <xsl:call-template name="quotify">
            <xsl:with-param name="string" select="$productName"/>
          </xsl:call-template>
          <xsl:text>','</xsl:text>
          <xsl:call-template name="quotify">
            <xsl:with-param name="string" select="label"/>
          </xsl:call-template>
          <xsl:text>','</xsl:text>
          <xsl:call-template name="quotify">
            <xsl:with-param name="string" select="file"/>
          </xsl:call-template>
          <xsl:text>')</xsl:text>
        </xsl:otherwise>
      </xsl:choose>
    </xsl:variable>

    <xsl:variable name="linkPath">
      <xsl:call-template name="checkLinkPath">
        <xsl:with-param name="linkPath" select="$link"/>
        <xsl:with-param name="rootPath" select="$demosroot"/>
      </xsl:call-template>
    </xsl:variable> 

    <li id="{$link}" class="node_element">
      <div class="node_element">
        <div class="wrapper img_wrapper"></div>
        <div class="wrapper text_wrapper">
          <a href="{$linkPath}" class="corrected_url"><xsl:value-of select="label"/></a>
        </div>
      </div>
    </li>
  </xsl:template>

  <xsl:template match="dependency">
    <xsl:if test="position()!=1">, </xsl:if>
    <xsl:value-of select="."/>
  </xsl:template>

  <xsl:template name="quotify">
    <xsl:param name="string"/>
    <xsl:variable name="percentEscaped">
      <xsl:call-template name="globalReplace">
        <xsl:with-param name="outputString" select="$string"/>
        <xsl:with-param name="target" select="'%'"/>
        <xsl:with-param name="replacement" select="'%25'"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="doubleUnQuoted">
      <xsl:call-template name="globalReplace">
        <xsl:with-param name="outputString" select="$percentEscaped"/>
        <xsl:with-param name="target" select="'&quot;'"/>
        <xsl:with-param name="replacement" select="'&amp;quot;'"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:call-template name="globalReplace">
      <xsl:with-param name="outputString" select="$doubleUnQuoted"/>
      <xsl:with-param name="target" select='"&apos;"'/>
      <xsl:with-param name="replacement" select='"&apos;&apos;"'/>
    </xsl:call-template>
  </xsl:template>

  <!-- Search and replace -->
  <!-- From http://www.xml.com/lpt/a/2002/06/05/transforming.html -->
  <xsl:template name="globalReplace">
    <xsl:param name="outputString"/>
    <xsl:param name="target"/>
    <xsl:param name="replacement"/>
    <xsl:choose>
      <xsl:when test="contains($outputString,$target)">
        <xsl:value-of select="concat(substring-before($outputString,$target),$replacement)"/>
        <xsl:call-template name="globalReplace">
          <xsl:with-param name="outputString" select="substring-after($outputString,$target)"/>
          <xsl:with-param name="target" select="$target"/>
          <xsl:with-param name="replacement" select="$replacement"/>
        </xsl:call-template>
      </xsl:when>
      <xsl:otherwise>
        <xsl:value-of select="$outputString"/>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

  <!-- A helper method used to build a localized path from a non-localized path.
       This method is called only when $languageDir is not '' and we know that a
       localized directory should exist. -->
  <xsl:template name="buildLocalizedPath">
    <xsl:param name="filepath"/>
    <xsl:choose>
      <xsl:when test="contains($filepath,'/')">
        <xsl:value-of select="concat(substring-before($filepath,'/'),'/')"/>
        <xsl:call-template name="buildLocalizedPath">
          <xsl:with-param name="filepath" select="substring-after($filepath,'/')"/>
        </xsl:call-template>
      </xsl:when>
      <xsl:otherwise>
        <xsl:value-of select="$languageDir"/>/<xsl:value-of select="$filepath"/>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

  <!-- A helper function that checks if a localized file path exists.  If it
       does, that file path will be returned, otherwise a specified default
       path will be. -->
  <xsl:template name="checkLocalizedFile">
    <xsl:param name="localizedFile"/>
    <xsl:param name="defaultFile"/>
    <xsl:choose>
      <xsl:when test="$localizedFile != ''">
        <xsl:variable name="fullpath" select="concat($demosroot,'/',$localizedFile)"/>
        <xsl:choose>
          <xsl:when test="file:exists(file:new($fullpath))">
            <xsl:value-of select="$localizedFile"/>
          </xsl:when>
          <xsl:otherwise>
            <xsl:value-of select="$defaultFile"/>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:when>
      <xsl:otherwise>
        <xsl:value-of select="$defaultFile"/>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

  <!-- A helper function to check the link path. If the link path is not a 
       matlab: link, pre-prend the root path. -->
  <xsl:template name="checkLinkPath">
    <xsl:param name="linkPath"/>    
    <xsl:param name="rootPath"/>
    <xsl:choose>
      <xsl:when test="contains($linkPath,'matlab:')"><xsl:value-of select="$linkPath"/></xsl:when>      
      <xsl:otherwise><xsl:value-of select="concat($rootPath,'/',$linkPath)"/></xsl:otherwise>      
    </xsl:choose>    
  </xsl:template>
  
</xsl:stylesheet>